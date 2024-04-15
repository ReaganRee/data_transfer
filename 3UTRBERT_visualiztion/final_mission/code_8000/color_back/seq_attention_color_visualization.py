import torch
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import numpy as np
from Bio import SeqIO
from tqdm import tqdm
np.set_printoptions(threshold=np.inf)

from transformers import BertModel, DNATokenizer
# from process_pretrain_data import get_kmer_sentence



def get_kmer_sentence(original_string, kmer=1, stride=1):
    # if kmer == -1:
    #     return original_string
    #
    # sentence = ""
    # original_string = original_string.replace("\n", "")
    # i = 0
    # while i < len(original_string)-kmer:
    #     sentence += original_string[i:i+kmer] + " "
    #     i += stride

    # return sentence[:].strip("\"")
    if kmer == -1:
        return original_string

    sequence = []
    original_string = original_string.replace("\n", "")
    for i in range(len(original_string) - kmer):
        sequence.append(original_string[i:i + kmer])

    sequence.append(original_string[-kmer:])
    return sequence

def format_attention(attention):
    squeezed = []
    for layer_attention in attention:
        # 1 x num_heads x seq_len x seq_len (each layer_attention)
        #print("layer_attention.shape: ", layer_attention.shape)
        if len(layer_attention.shape) != 4:
            raise ValueError("The attention tensor does not have the correct number of dimensions. Make sure you set "
                             "output_attentions=True when initializing your model.")
        squeezed.append(layer_attention.squeeze(0)) # layer_attention.squeeze(0), remove the 0th dimension if it is 1, num_heads x seq_len x seq_len
        #print("squeeze:: ", layer_attention.squeeze(0).shape)
    # num_layers x num_heads x seq_len x seq_len
    #print(torch.stack(squeezed).shape)
    return torch.stack(squeezed) #combine all metrics in squeezed to one 12 x 12 x seq_len x seq_len

def get_attention_dna(model, tokenizer, sentence_a, start, end):
    inputs = tokenizer.encode_plus(sentence_a, sentence_b=None, return_tensors='pt', add_special_tokens=True, max_length=512)
    input_ids = inputs['input_ids']
    #print("input_ids: ", input_ids)
    attention = model(input_ids)[-1]
    #print("attention len: ", len(attention))
    #print("attention shape: ", attention[-1].shape)
    input_id_list = input_ids[0].tolist() # Batch index 0
    #print("input_id_list: ", input_id_list)
    tokens = tokenizer.convert_ids_to_tokens(input_id_list)
    print("tokens: ", tokens)
    attn = format_attention(attention)
    # print(attn)# num_layers x num_heads x seq_len x seq_len
    attn_score = []
    for i in range(1, len(tokens)-1):
        help = attn[start:end+1,:,0,i]
        # print(help)
        attn_score.append(float(attn[start:end+1,:,0,i].sum())) # 此处start:end+1限制取atten的范围在11 layer,逻辑是每个head中CLS对位置i的token的和(reconfirm this!!!)
    return attn_score

def get_real_score(attention_scores, kmer, metric):
    counts = np.zeros([len(attention_scores)+kmer-1])
    real_scores = np.zeros([len(attention_scores)+kmer-1])

    if metric == "mean":
        for i, score in enumerate(attention_scores):
            for j in range(kmer):
                counts[i+j] += 1.0
                real_scores[i+j] += score

        real_scores = real_scores/counts
    else:
        pass

    return real_scores




SEQUENCE = ""

def generate_attention_average(args):

    if args.kmer == 0:
        KMER_LIST = [3] # [3,4,5,6]

        for kmer in KMER_LIST:
            tokenizer_name = 'dna' + str(kmer)
            model_path = os.path.join(args.model_path, str(kmer))
            model = BertModel.from_pretrained(model_path, output_attentions=True)
            tokenizer = DNATokenizer.from_pretrained(tokenizer_name, do_lower_case=False)
            raw_sentence = args.sequence if args.sequence else SEQUENCE
            sentence_a = ' '.join(get_kmer_sentence(raw_sentence, kmer))
            tokens = sentence_a.split()

            attention = get_attention_dna(model, tokenizer, sentence_a, start=args.start_layer, end=args.end_layer)
            attention_scores = np.array(attention).reshape(np.array(attention).shape[0],1)
            # attention_scores[0] = 0

            real_scores = get_real_score(attention_scores, kmer, args.metric)
            real_scores = real_scores / np.linalg.norm(real_scores)

            if kmer != KMER_LIST[0]:
                scores += real_scores.reshape(1, real_scores.shape[0])
            else:
                scores = real_scores.reshape(1, real_scores.shape[0])

    else:
        # load model and calculate attention
        #tokenizer_name = 'dna' + str(args.kmer)
        model_path = args.model_path
        model = BertModel.from_pretrained(model_path, output_attentions=True)
        tokenizer = DNATokenizer.from_pretrained("/Users/reagan/Desktop/3UTRBERT_visualiztion/saliency_based_color_background/3-new-12w-0_a/", do_lower_case=False)
        #print(tokenizer)
        raw_sentence = args.sequence if args.sequence else SEQUENCE
        sentence_a = ' '.join(get_kmer_sentence(raw_sentence, args.kmer))
        #print(sentence_a)
        tokens = sentence_a.split()
        attention = get_attention_dna(model, tokenizer, sentence_a, start=args.start_layer, end=args.end_layer)
        #print(np.array(attention).shape)
        #print("reagan want to see: ", np.array(attention).shape[0])
        attention_scores = np.array(attention).reshape(np.array(attention).shape[0],1)
        #print("reagan want to see: ", attention_scores.shape)

        # attention_scores[0] = 0

        real_scores = get_real_score(attention_scores, args.kmer, args.metric)
        scores = real_scores.reshape(1, real_scores.shape[0])

    ave = np.sum(scores)/scores.shape[1]
    # print(ave)
    # print(scores)
    return real_scores


def highlighter(i, weight, input_seq_list):
    help = weight[i]
    color = '#%02X%02X%02X' % (
        255, int(255*(1 - weight[i])), int(255*(1 - weight[i])))
    word = '<span style="background-color:' +color+ '">' +input_seq_list[i]+ '</span>'
    return word





def main(sequence, model_path):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kmer",
        default=3,
        type=int,
        help="K-mer",
    )
    parser.add_argument(
        "--model_path",
        default=model_path,
        type=str,
        help="The path of the finetuned model",
    )
    parser.add_argument(
        "--start_layer",
        default=11,
        type=int,
        help="Which layer to start",
    )
    parser.add_argument(
        "--end_layer",
        default=11,
        type=int,
        help="which layer to end",
    )
    parser.add_argument(
        "--metric",
        default="mean",
        type=str,
        help="the metric used for integrate predicted kmer result to real result",
    )
    parser.add_argument(
        "--sequence",
        default=sequence,
        type=str,
        help="the sequence for visualize",
    )

    args = parser.parse_args()
    weight = generate_attention_average(args)
    # print("average attention: ", weight)
    # print(type(weight))
    # weight[0] = 0.9
    # weight[1] = 0.8
    # weight[2] = 0.61
    # weight[6] = 0.1
    # weight[7] = 0.23
    # weight[8] = 0.43
    # weight[9] = 0.3
    # weight[10] = 0.26
    # weight[11] = 0.18
    # weight[12] = 0.08

    text = ''.join([highlighter(i, weight, sequence) for i in range(len(sequence))])

    #display(HTML(text))
    # with open("/Users/reagan/Desktop/3UTRBERT_visualiztion/colored_background_visualization/visualization.html", "w") as file:
    #     file.write(text)
    return weight.tolist()

if __name__ == "__main__":
    # attentions_dict = {}
    # records_list = list(SeqIO.parse(
    #     "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6a_3UTR_ind_pos.fa",
    #     "fasta"))
    # for seq_record in records_list[1:2]:
    #     length = len(str(seq_record.seq))
    #     #print(length)
    #     if length > 510:
    #         counter = length // 510
    #         start = 0
    #         catcher = []
    #         while counter > 0:
    #             sequence = str(seq_record.seq)[start:start+510]
    #             attention = main(sequence)
    #             catcher.extend(attention)
    #             start += 510
    #             counter -= 1
    #         catcher.extend(main(str(seq_record.seq)[start:]))
    #         print("cater_len: ", len(catcher))
    #         attentions_dict[seq_record.id] = catcher
    #
    #     elif length > 5 and len(str(seq_record.seq)) <= 510:
    #         attention2 = main(str(seq_record.seq))
    #         attentions_dict[seq_record.id] = attention2
    # np.save("/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/m6a_attentions_1_for_test", attentions_dict)
    for cell_line in ["A549", "CD8T", "ESC", "HCT116", "HEK293", "HEK293T", "Hela", "HepG2", "MOLM13"]:
        model_path = "/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/m6a_9/{}".format(cell_line)
        with open("/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/{}.txt".format(cell_line), "r") as file:
            seq_list = file.readlines()
            for seq in seq_list:
                sequence = seq.strip()
                attention = main(sequence, model_path)


    sequence = "GCUGGAUAACUUAUUUAUGGACUGUUGGGGAUGAGAGCAGG"
    attention = main(sequence, model_path)
    #print(attention)


    # result = []
    # i = 0
    # for seq_record in tqdm(SeqIO.parse(
    #     "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/for_test.fa",
    #     "fasta")):
    #     # print(seq_record.seq)
    #     i += 1
    #     print(i)
        # if i > 0:
        #     result.append(seq_record)
        #     i -= 1
        # else:
        #     SeqIO.write(result, "/Users/reagan/Desktop/Research/UT_Zhaolei_Zhang_Lab/Bert_visulization/visualization_4/for_test.fa", 'fasta')

