import os, sys

import numpy
import numpy as np
import matplotlib as mpl
from PIL import Image
mpl.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from scipy.misc import imresize

package_directory = os.path.dirname(os.path.abspath(__file__))
#print("package_directory: ", package_directory)
acgu_path = os.path.join(package_directory, 'acgu.npz')
chars = np.load("/Users/reagan/Desktop/3UTRBERT_visualiztion/saliency_based_color_background/acgu.npz", allow_pickle=True)['data']


def normalize_pwm(pwm, factor=None, MAX=None):
    if MAX is None:
        MAX = np.max(np.abs(pwm))
    pwm = pwm / MAX
    if factor:
        pwm = np.exp(pwm * factor)
    norm = np.outer(np.ones(pwm.shape[0]), np.sum(np.abs(pwm), axis=0))
    return pwm / norm


def get_nt_height(pwm, height, norm):

    def entropy(p):
        s = 0
        for i in range(len(p)):
            if p[i] > 0:
                s -= p[i] * np.log2(p[i])
        return s

    num_nt, num_seq = pwm.shape
    heights = np.zeros((num_nt, num_seq))
    for i in range(num_seq):
        if norm == 1:
            total_height = height
        else:
            total_height = (np.log2(num_nt) - entropy(pwm[:, i])) * height

        heights[:, i] = np.floor(
            pwm[:, i] * np.minimum(total_height, height * 2))

    return heights.astype(int)


def seq_logo(pwm, height=30, nt_width=10, norm=0, alphabet='rna',
             colormap='standard'):

    heights = get_nt_height(pwm, height, norm)
    num_nt, num_seq = pwm.shape
    width = np.ceil(nt_width * num_seq).astype(int)

    max_height = height*2
    logo = np.ones((max_height, width, 3)).astype(int) * 255
    for i in range(num_seq):
        nt_height = np.sort(heights[:, i])
        index = np.argsort(heights[:, i])
        remaining_height = np.sum(heights[:, i])
        offset = max_height - remaining_height

        for j in range(num_nt):
            if nt_height[j] <= 0:
                continue
            # resized dimensions of image

            nt_img = numpy.array(Image.fromarray(chars[index[j]]).resize(size = (nt_width, nt_height[j])))#imresize(chars[index[j]], (nt_height[j], nt_width))
            # determine location of image
            height_range = range(remaining_height - nt_height[j],
                                 remaining_height)
            width_range = range(i * nt_width, i * nt_width + nt_width)
            # 'annoying' way to broadcast resized nucleotide image
            if height_range:
                for k in range(3):
                    for m in range(len(width_range)):
                        logo[height_range + offset, width_range[m], k] = nt_img[
                                                                         :, m,
                                                                         k]

            remaining_height -= nt_height[j]

    return logo.astype(np.uint8)


def plot_saliency(X, W, outdir, nt_width=100, norm_factor=3, str_null=None
                  ):
    # filter out zero-padding
    plot_index = np.where(np.sum(X[:4, :], axis=0) != 0)[0]
    #print("plot_index: ", plot_index)
    num_nt = len(plot_index)
    trace_width = num_nt * nt_width
    trace_height = 400

    seq_str_mode = False
    if X.shape[0] > 4:
        seq_str_mode = True
        assert str_null is not None, "Null region is not provided."

    # sequence logo
    img_seq_raw = seq_logo(X[:4, plot_index], height=nt_width,
                           nt_width=nt_width)

    if seq_str_mode:
        # structure line
        str_raw = X[4, plot_index]
        if str_null.sum() > 0:
            str_raw[str_null.T == 1] = -0.01

        line_str_raw = np.zeros(trace_width)
        for v in range(str_raw.shape[0]):
            line_str_raw[v * nt_width:(v + 1) * nt_width] = (1 - str_raw[
                v]) * trace_height
            # i+=1

    # sequence saliency logo
    seq_sal = normalize_pwm(W[:4, plot_index], factor=norm_factor)
    img_seq_sal_logo = seq_logo(seq_sal, height=nt_width * 5, nt_width=nt_width)
    img_seq_sal = numpy.array(Image.fromarray(W[:4, plot_index]).resize((trace_width, trace_height)))#imresize(W[:4, plot_index], size=(trace_height, trace_width))

    if seq_str_mode:
        # structure saliency logo
        str_sal = W[4, plot_index].reshape(1, -1)
        img_str_sal = numpy.array(Image.fromarray(str_sal).resize((trace_width, trace_height)))#imresize(str_sal, size=(trace_height, trace_width))

    # plot
    fig = plt.figure(figsize=(10.1, 2))
    gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[2.5, 1, 0.5, 1])
    cmap_reversed = mpl.cm.get_cmap('jet')

    ax = fig.add_subplot(gs[0, 0])
    ax.axis('off')
    ax.imshow(img_seq_sal_logo)
    plt.text(x=trace_width - 400, y=10, s='PrismNet', fontsize=4)

    ax = fig.add_subplot(gs[1, 0])
    ax.axis('off')
    ax.imshow(img_seq_sal, cmap=cmap_reversed)

    ax = fig.add_subplot(gs[2, 0])
    ax.axis('off')
    ax.imshow(img_seq_raw)

    if seq_str_mode:
        ax = fig.add_subplot(gs[3, 0])
        ax.axis('off')
        ax.imshow(img_str_sal, cmap=cmap_reversed)
        ax.plot(line_str_raw, '-', color='r', linewidth=1, scalex=False,
                scaley=False)

        # plot balck line to hide the -1(NULL structure score)
        x = (np.zeros(trace_width) + (1 + 0.01)) * trace_height + 1.5
        ax.plot(x, '-', color='white', linewidth=1.2, scalex=False,
                scaley=False)

    plt.subplots_adjust(wspace=0, hspace=0)

    # save figure
    filepath = outdir
    fig.savefig(filepath, format='pdf', dpi=2000, bbox_inches='tight')
    plt.close('all')

if __name__ == "__main__":
    from color_back.seq_attention_color_visualization import main
    for cell_line in ["A549", "CD8T", "ESC", "HCT116", "HEK293", "HEK293T", "Hela", "HepG2", "MOLM13"]:
        print(cell_line)
        model_path = "/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/m6a_9/{}".format(cell_line)
        with open("/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/sequences/{}.txt".format(cell_line), "r") as file:
            seq_list = file.readlines()
            for seq in seq_list[:20]:
                sequence1 = seq.strip()
                attention1 = main(sequence1, model_path)
                attention = attention1[2:] + attention1[:2]
                X1 = []
                W1 = []
                sequence = list(sequence1)
                #print(type(X1))
                for i in range(len(sequence)):
                    sub_list = [0, 0, 0, 0]
                    sub_list_2 = [0, 0, 0, 0]
                    if sequence[i] == "A":
                        sub_list[0] = 1
                        sub_list_2[0] = attention[i]

                    if sequence[i] == "C":
                        sub_list[1] = 1
                        sub_list_2[1] = attention[i]

                    if sequence[i] == "G":
                        sub_list[2] = 1
                        sub_list_2[2] = attention[i]

                    if sequence[i] == "U":
                        sub_list[3] = 1
                        sub_list_2[3] = attention[i]
                    #print(type(X1))
                    X1.append(sub_list)
                    W1.append(sub_list_2)
                    X2 = np.array(X1)

                    W2 = np.array(W1)

                plot_saliency(X2.T, W2.T, "/Users/reagan/Desktop/final_mission/saliency_m6a_9cellines/result_2000_pdf/{cell}_{seq}.pdf".format(cell=cell_line,seq=sequence1))






    # old_W1 =[[0, 0, 0, 0.4541873335838318], [0, 0, 0, 0.25428760424256325], [0, 0, 0, 0.18189162760972977], [0, 0, 0, 0.032440901113053165], [0, 0, 0.018417546525597572, 0], [0, 0, 0.01160734953979651, 0], [0, 0.01596746624757846, 0, 0], [0, 0.018139773979783058, 0, 0], [0, 0, 0.026083222900827725, 0], [0.0238228611027201, 0, 0, 0], [0, 0, 0.02496358100324869, 0], [0, 0, 0.023182273842394352, 0], [0, 0, 0.026698917771379154, 0], [0.026187484463055927, 0, 0, 0], [0.01777716384579738, 0, 0, 0], [0, 0, 0.021110464197893936, 0], [0, 0.023004585566620033, 0, 0], [0.02281481431176265, 0, 0, 0], [0, 0, 0.6003079107031226, 0], [0, 0, 1.2627105865006645, 0], [1.4850315650304158, 0, 0, 0], [0, 1.903178056081136, 0, 0], [1.7600653966267903, 0, 0, 0], [0, 0, 1.8821558554967244, 0], [0.8820997371027867, 0, 0, 0], [0, 0, 0.37661008226374787, 0], [0.041938296829660736, 0, 0, 0], [0.03837632077435652, 0, 0, 0], [0, 0, 0.02349892444908619, 0], [0, 0, 0.01639996903638045, 0], [0.017516742770870525, 0, 0, 0], [0, 0, 0.019240723301966984, 0], [0, 0.03326437125603358, 0, 0], [0.03172354958951473, 0, 0, 0], [0, 0, 0.10306435823440552, 0], [0, 0, 0.09380561485886574, 0], [0.09976128054161866, 0, 0, 0], [0.02195730184515317, 0, 0, 0], [0, 0, 0.02162607138355573, 0], [0, 0.01976731512695551, 0, 0], [0, 0, 0, 0.03057834506034851]]


    #
    # for j in range(len(old_W1)):
    #     for i in range(4):
    #         if j == 18 or j==19:
    #             old_W1[j][i] = old_W1[j][i] * 2


# W1[0] = 0.9
# W1[1] = 0.8
# W1[2] = 0.61
# W1[6] = 0.1
# W1[7] = 0.23
# W1[8] = 0.43
# W1[9] = 0.3
# W1[10] = 0.26
# W1[11] = 0.18
# W1[12] = 0.08
