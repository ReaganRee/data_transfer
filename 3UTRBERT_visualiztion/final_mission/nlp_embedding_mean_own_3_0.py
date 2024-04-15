import os
os.environ['CUDA_VISIBLE_DEVICES']='0'

import torch
import numpy as np
import pandas as pd

from Bio import SeqIO
from torch import cuda
from torch.utils.data import DataLoader, Dataset
from transformers import BertTokenizer, BertModel, BertConfig
from keras.preprocessing.sequence import pad_sequences

#from color_back.seq_attention_color_visualization import main

def mk_dir(dir):
    try:
        os.makedirs(dir)
    except OSError:
        print('Can not make directory:', dir)

class ChunkDataset(Dataset):
    def __init__(self, text, labels, tokenizer, chunk_len=512, overlap_len=0):
        self.tokenizer = tokenizer
        self.text = text
        self.labels = labels
        self.overlap_len = overlap_len
        self.chunk_len = chunk_len

    def __len__(self):
        return len(self.labels)

    def chunk_tokenizer(self, tokenized_data, targets):
        input_ids_list = []
        attention_mask_list = []
        token_type_ids_list = []
        targets_list = []

        previous_input_ids = tokenized_data["input_ids"]
        previous_attention_mask = tokenized_data["attention_mask"]
        previous_token_type_ids = tokenized_data["token_type_ids"]
        remain = tokenized_data.get("overflowing_tokens")


        input_ids_list.append(torch.tensor(previous_input_ids, dtype=torch.long))
        attention_mask_list.append(torch.tensor(previous_attention_mask, dtype=torch.long))
        token_type_ids_list.append(torch.tensor(previous_token_type_ids, dtype=torch.long))
        targets_list.append(torch.tensor(targets, dtype=torch.long))


        if remain:  # if there is any overflowing tokens
            # remain = torch.tensor(remain, dtype=torch.long)
            idxs = range(len(remain) + self.chunk_len)
            idxs = idxs[(self.chunk_len - self.overlap_len - 2)
                        ::(self.chunk_len - self.overlap_len - 2)]
            input_ids_first_overlap = previous_input_ids[-(self.overlap_len + 1):-1]

            start_token = [1]
            end_token = [2]

            for i, idx in enumerate(idxs):
                if i == 0:
                    input_ids = input_ids_first_overlap + remain[:idx]
                elif i == len(idxs):
                    input_ids = remain[idx:]
                elif previous_idx >= len(remain):
                    break
                else:
                    input_ids = remain[(previous_idx - self.overlap_len):idx]

                previous_idx = idx

                nb_token = len(input_ids) + 2
                attention_mask = np.ones(self.chunk_len)
                attention_mask[nb_token:self.chunk_len] = 0
                token_type_ids = np.zeros(self.chunk_len)
                input_ids = start_token + input_ids + end_token
                if self.chunk_len - nb_token > 0:
                    padding = np.zeros(self.chunk_len - nb_token)
                    input_ids = np.concatenate([input_ids, padding])

                input_ids_list.append(torch.tensor(input_ids, dtype=torch.long))
                attention_mask_list.append(torch.tensor(attention_mask, dtype=torch.long))
                token_type_ids_list.append(torch.tensor(token_type_ids, dtype=torch.long))
                targets_list.append(torch.tensor(targets, dtype=torch.long))

        return ({
            'ids': input_ids_list,
            'mask': attention_mask_list,
            'token_type_ids': token_type_ids_list,
            'targets': targets_list,
            'len': [torch.tensor(len(targets_list), dtype=torch.long)]
        })

    def __getitem__(self, index):
        text = " ".join(str(self.text[index]).split())
        targets = self.labels[index]

        data = self.tokenizer.encode_plus(
            text=text,
            text_pair=None,
            add_special_tokens=True,
            max_length=self.chunk_len,
            truncation=True,
            pad_to_max_length=True,
            return_token_type_ids=True,
            return_overflowing_tokens=True
        )

        chunk_token = self.chunk_tokenizer(data, targets)
        return chunk_token

def chunk_collate_fn(batches):
    """
    Create batches for ChunkDataset
    """
    return [{key: torch.stack(value) for key, value in batch.items()} for batch in batches]


class MyDataset(Dataset):
    def __init__(self,df):
        self.X = df['SequenceID'].to_list()
        self.Y = df['Label']
    def __len__(self):
        return len(self.X)
    def __getitem__(self,index):
        return self.X[index], self.Y.iloc[index]

def load_seq(seq_path, kmer):
    sequence = []
    seq_label = []
    seq_content_list = []

    for seq_record in SeqIO.parse(seq_path, "fasta"):
        seq_label.append(str((seq_record.id).split(',')[0]))
        seq_origin = str(seq_record.seq.strip())
        seq_origin = seq_origin.upper().replace('T', 'U')
        sequence.append(seq2kmer(seq_origin.strip(), kmer))
    df = pd.DataFrame(data={'SequenceID': sequence, 'Label': seq_label})
    return df

def seq2kmer(seq, k):
    """
    Convert original sequence to kmers

    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.

    Returns:
    kmers -- str, kmers separated by space
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    #kmer = re.findall(r'\w{3}',seq)
    kmers = " ".join(kmer)
    return kmers

def kmer2seq(kmers):
    """
    Convert kmers to original sequence

    Arguments:
    kmers -- str, kmers separated by space.

    Returns:
    seq -- str, original sequence.

    """
    kmers_list = kmers.split(" ")
    bases = [kmer[0] for kmer in kmers_list[0:-1]]
    bases.append(kmers_list[-1])
    seq = "".join(bases)
    assert len(seq) == len(kmers_list) + len(kmers_list[0]) - 1
    return seq

def vectorize_labels(all_labels):
    """
    :return: dict of vectorized labels per split and total number of labels
    """
    result = {}
    for split in all_labels:
        result[split] = np.array(all_labels[split])
    return result

def remove_special_token(embedding, attention_mask):
    transform = []
    for seq_num in range(len(embedding)):
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][1:seq_len-1]
        transform.append(seq_emd)
    transform_emb = np.vstack(transform)
    return transform_emb

def prepare_data(data_path, dataset_num, num_labels, kmer):
    """
    return: dicts of lists of documents and labels and number of labels
    """
    if not os.path.exists(data_path):
        raise Exception("Data path not found: {}".format(data_path))

    text_set = {'Train_fold': [], 'Validation_fold': [], 'Test_fold': []}
    label_set = {'Train_fold': [], 'Validation_fold': [], 'Test_fold': []}

    dataset_split = ['Train_fold', 'Validation_fold', 'Test_fold']
    for each_item in dataset_split:
        wholepath = data_path + each_item + str(dataset_num) + '.fasta'
        df_dataset = load_seq(wholepath, kmer)
        df_dataset["Label"] = df_dataset["Label"].apply(lambda x: list(map(int, x)))
        df_dataset = MyDataset(df_dataset)
        for item in df_dataset:
            #print(item[0])
            text_set[each_item].append(item[0])
            label_set[each_item].append(item[1])

    vectorized_labels = vectorize_labels(label_set)
    return text_set, vectorized_labels, num_labels

def create_dataloader(dataset_class, text_set, label_set, tokenizer, max_length, batch_size, num_workers):
    """
    Create appropriate dataloaders for the given data
    """
    dataloaders = {}

    if 'Train_fold' in text_set.keys():
        split = 'Train_fold'
        dataset = dataset_class(text_set[split], label_set[split], tokenizer, max_length)
        if isinstance(dataset, ChunkDataset):
            dataloaders[split] = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers,
                                            pin_memory=True, collate_fn=chunk_collate_fn)
        else:
            dataloaders[split] = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)

    for split in ['Validation_fold', 'Test_fold']:
        dataset = dataset_class(text_set[split], label_set[split], tokenizer, max_length)
        if isinstance(dataset, ChunkDataset):
            dataloaders[split] = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers,
                                            pin_memory=True, collate_fn=chunk_collate_fn)
        else:
            dataloaders[split] = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, pin_memory=True)

    return dataloaders

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


if __name__ == "__main__":
    data_path = '/code_8000/5_fold_data_final/'
    output_path = '/code_8000/output/'
    dataset_num= 0
    classes = 7
    kmer = 3
    max_length = 512
    batch_size = 5
    num_workers = 0
    fixed_length = 8000

    mk_dir(output_path)

    tokenizer = BertTokenizer.from_pretrained('/Users/reagan/Desktop/final_mission/code_8000/3-new-12w-0/', do_lower_case=False)
    text_set, label_set, num_labels = prepare_data(data_path, dataset_num, classes, kmer)

    dataset_class = ChunkDataset
    dataloaders = create_dataloader(dataset_class, text_set, label_set, tokenizer, max_length, batch_size, num_workers)


    model = BertModel.from_pretrained('/code_8000/3-new-12w-0/', config=BertConfig.from_pretrained('/Users/reagan/Desktop/final_mission/code_8000/3-new-12w-0/', output_attentions=True))

    #model = BertModel.from_pretrained("/home/wangyansong/mRNA/3-new-12w-0/",output_hidden_states=True)
    device = 'cuda' if cuda.is_available() else 'cpu'
    model = model.to(device)
    model = model.eval()

    train_attn = []
    valid_attn = []
    test_attn = []


    for each_id in  dataloaders.keys():
        if each_id == 'Train_fold':
            split = 'Train_fold'
            with torch.no_grad():
                for batch_idx, data in enumerate(dataloaders[split], 0):
                    for each_item in data:
                        ids = each_item['ids'].to(device, dtype = torch.long)
                        print(ids.shape)
                        #mask = each_item['mask'].to(device, dtype = torch.long)
                        token_type_ids = each_item['token_type_ids'].to(device, dtype = torch.long)
                        outputs = model(input_ids=ids)
                        attention = outputs[-1]
                        # print("attention len: ", len(attention))
                        # print("attention shape: ", attention[-1].shape)
                        input_id_list = ids[0].tolist()  # Batch index 0
                        # print("input_id_list: ", input_id_list)
                        tokens = tokenizer.convert_ids_to_tokens(input_id_list)
                        print("tokens: ", tokens)
                        attn = format_attention(attention)
                        # print(attn)# num_layers x num_heads x seq_len x seq_len
                        attn_score = []
                        for i in range(1, len(tokens) - 1):
                            help = attn[11:11 + 1, :, 0, i]
                            # print(help)
                            attn_score.append(
                                float(attn[11:11 + 1, :, 0, i].sum()))
                            # hidden_states = outputs[2]
                        attention_2 = attn_score
                        attention_scores = np.array(attention_2).reshape(
                            np.array(attention_2).shape[0], 1)
                        train_real_scores = get_real_score(attention_scores,
                                                     3, "mean")
                        print("train_shape: ", train_real_scores.shape)
                        train_attn.extend(train_real_scores)
                        print("train_attn: ", len(train_attn))

                        # hidden_states = (hidden_states[-1] + hidden_states[1]).cpu().numpy()
                        # transform_emb = remove_special_token(hidden_states, mask)
                        # embedding_pad = np.pad(transform_emb, ((0,fixed_length-transform_emb.shape[0]),(0,0)), 'mean')
                        # embedding_pad = np.mean(embedding_pad, axis=1)
                        # train_embedding.append(embedding_pad[:, np.newaxis])

        elif each_id == 'Validation_fold':
            split = 'Validation_fold'
            with torch.no_grad():
                for batch_idx, data in enumerate(dataloaders[split], 0):
                    for each_item in data:
                        ids = each_item['ids'].to(device, dtype = torch.long)
                        #mask = each_item['mask'].to(device, dtype = torch.long)
                        token_type_ids = each_item['token_type_ids'].to(device, dtype = torch.long)
                        outputs = model(input_ids=ids)
                        attention = outputs[-1]
                        # print("attention len: ", len(attention))
                        # print("attention shape: ", attention[-1].shape)
                        input_id_list = ids[0].tolist()  # Batch index 0
                        # print("input_id_list: ", input_id_list)
                        tokens = tokenizer.convert_ids_to_tokens(input_id_list)
                        #print("tokens: ", tokens)
                        attn = format_attention(attention)
                        # print(attn)# num_layers x num_heads x seq_len x seq_len
                        attn_score = []
                        for i in range(1, len(tokens) - 1):
                            help = attn[11:11 + 1, :, 0, i]
                            # print(help)
                            attn_score.append(
                                float(attn[11:11 + 1, :, 0, i].sum()))
                            # hidden_states = outputs[2]
                        attention_2 = attn_score
                        attention_scores = np.array(attention_2).reshape(
                            np.array(attention_2).shape[0], 1)
                        validation_real_scores = get_real_score(attention_scores,
                                                     3, "mean")
                        print("vali_shape: ", validation_real_scores.shape)
                        valid_attn.extend(validation_real_scores)
                        print("vali_attn: ", len(valid_attn))
                        #scores = real_scores.reshape(1, real_scores.shape[0])
                        # hidden_states = outputs[2]
                        # ##hidden_states = model(input_ids=ids,attention_mask=mask, output_hidden_states=True).hidden_states
                        # #outputs = model(input_ids=ids,attention_mask=mask)
                        # #print(outputs)
                        # hidden_states = (hidden_states[-1] + hidden_states[1]).cpu().numpy()
                        # transform_emb = remove_special_token(hidden_states, mask)
                        # embedding_pad = np.pad(transform_emb, ((0,fixed_length-transform_emb.shape[0]),(0,0)), 'mean')
                        # embedding_pad = np.mean(embedding_pad, axis=1)
                        # valid_embedding.append(embedding_pad[:, np.newaxis])
        else:
            split = 'Test_fold'
            with torch.no_grad():
                for batch_idx, data in enumerate(dataloaders[split], 0):
                    for each_item in data:
                        ids = each_item['ids'].to(device, dtype = torch.long)
                        #mask = each_item['mask'].to(device, dtype = torch.long)
                        token_type_ids = each_item['token_type_ids'].to(device, dtype = torch.long)
                        outputs = model(input_ids=ids)
                        attention = outputs[-1]
                        # print("attention len: ", len(attention))
                        # print("attention shape: ", attention[-1].shape)
                        input_id_list = ids[0].tolist()  # Batch index 0
                        # print("input_id_list: ", input_id_list)
                        tokens = tokenizer.convert_ids_to_tokens(input_id_list)
                        #print("tokens: ", tokens)
                        attn = format_attention(attention)
                        # print(attn)# num_layers x num_heads x seq_len x seq_len
                        attn_score = []
                        for i in range(1, len(tokens) - 1):
                            help = attn[11:11 + 1, :, 0, i]
                            # print(help)
                            attn_score.append(
                                float(attn[11:11 + 1, :, 0, i].sum()))
                            # hidden_states = outputs[2]
                        attention_2 = attn_score
                        attention_scores = np.array(attention_2).reshape(
                            np.array(attention_2).shape[0], 1)
                        test_real_scores = get_real_score(attention_scores,
                                                     3, "mean")
                        print("test_shape: ", test_real_scores.shape)
                        test_attn.extend(test_real_scores)
                        print("test_attn: ", len(test_attn))
                        #scores = real_scores.reshape(1, real_scores.shape[0])
                        # hidden_states = outputs[2]
                        # #hidden_states = model(input_ids=ids,attention_mask=mask, return_dict=True, output_hidden_states=True).hidden_states
                        # hidden_states = (hidden_states[-1] + hidden_states[1]).cpu().numpy()
                        # transform_emb = remove_special_token(hidden_states, mask)
                        # embedding_pad = np.pad(transform_emb, ((0,fixed_length-transform_emb.shape[0]),(0,0)), 'mean')
                        # embedding_pad = np.mean(embedding_pad, axis=1)
                        # test_embedding.append(embedding_pad[:, np.newaxis])

    #print(np.array(train_embedding).shape)
    #print(np.array(valid_embedding).shape)
    #print(np.array(test_embedding).shape)
    #########################################################################
    #split = 'Train_fold'
    #for each_item in text_set[split]:
        #seq_length = len(kmer2seq(each_item))
        #print(seq_length)
        #train_length.append(seq_length)

    #for i1, i2 in zip(train_embedding, train_length):
        #delete_row = i1.shape[0]-i2
        #arr_process = np.delete(i1, np.s_[-delete_row: ], axis=0)
        #print(arr_process.shape)
    #########################################################################

    np.save(output_path + 'Train_fold' + str(dataset_num) + '.npy', np.array(train_attn))
    print("when saving, train_attn: ", len(train_attn))
    np.save(output_path + 'Validation_fold' + str(dataset_num) + '.npy', np.array(valid_attn))
    print("when saving, vali_attn: ", len(valid_attn))
    np.save(output_path + 'Test_fold' + str(dataset_num) + '.npy', np.array(test_attn))
    print("when saving, test_attn: ", len(test_attn))
