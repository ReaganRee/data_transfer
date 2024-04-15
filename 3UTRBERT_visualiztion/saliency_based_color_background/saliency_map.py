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
print("package_directory: ", package_directory)
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


def plot_saliency(X, W, nt_width=100, norm_factor=3, str_null=None,
                  outdir="/Users/reagan/Desktop/3UTRBERT_visualiztion/saliency_based_color_background/test_2.pdf"):
    # filter out zero-padding
    plot_index = np.where(np.sum(X[:4, :], axis=0) != 0)[0]
    print("plot_index: ", plot_index)
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
    fig.savefig(filepath, format='pdf', dpi=300, bbox_inches='tight')
    plt.close('all')

if __name__ == "__main__":
    old_W1 = [[0, 0, 0.4, 0], [0, 0.37, 0, 0], [0, 0, 0.305, 0], [0, 0, 0.1165863648056984, 0], [0, 0, 0, 0.1262977346777916], [0, 0, 0.1158574769894282, 0], [0, 0.1, 0, 0], [0, 0, 0.23, 0], [0.43, 0, 0, 0], [0, 0, 0.3, 0], [0, 0, 0.26, 0], [0, 0.18, 0, 0], [0, 0.08, 0, 0], [0, 0.16988486051559448, 0, 0], [0, 0, 0.17785650491714478, 0], [0, 0, 0.19783757130304971, 0], [0, 0.16283024350802103, 0, 0], [0, 0, 0.1610820988814036, 0], [0, 0, 0.12792385121186575, 0], [0, 0.16131782035032907, 0, 0], [0, 0, 0.16823380688826242, 0], [0, 0, 0, 0.1765719602505366], [0, 0.13361006478468576, 0, 0], [0, 0, 0.1288908471663793, 0], [0, 0, 0.19307716687520346, 0], [0, 0.2348264455795288, 0, 0], [0, 0, 0, 0.24266854425271353], [0, 0, 0.20591233174006143, 0], [0.1848839819431305, 0, 0, 0], [0, 0, 0.14753215511639914, 0], [0, 0, 0.11802399158477783, 0], [0, 0, 0.10346572597821553, 0], [0.10981037716070811, 0, 0, 0], [0, 0.11407755066951115, 0, 0], [0, 0.1546043778459231, 0, 0], [0.1585493658979734, 0, 0, 0], [0, 0, 0.18343166013558707, 0], [0, 0.12861687938372293, 0, 0], [0, 0, 0.14801567792892456, 0], [0, 0, 0.11622248589992523, 0], [0, 0, 0.16957386334737143, 0], [0.16878563662370047, 0, 0, 0], [0, 0.21812125047047934, 0, 0], [0, 0, 0, 0.21347567439079285], [0, 0, 0.2587843934694926, 0], [0, 0, 0.22777721285820007, 0], [0, 0, 0.18342783053716025, 0], [0, 0, 0.14753641684850058, 0], [0, 0.1256199354926745, 0, 0], [0, 0, 0.10477624957760175, 0], [0.06022505337993304, 0, 0, 0], [0.07112373039126396, 0, 0, 0], [0, 0.07341627031564713, 0, 0], [0, 0.0864507829149564, 0, 0], [0, 0.08575241019328435, 0, 0], [0, 0, 0.08578680704037349, 0], [0, 0, 0.09836067507664363, 0], [0, 0.08756613731384277, 0, 0], [0, 0, 0.102688267827034, 0], [0, 0, 0, 0.07496338834365208], [0, 0, 0.07266324013471603, 0], [0, 0, 0.059536270797252655, 0], [0, 0, 0.07699936131636302, 0], [0, 0.08458556234836578, 0, 0], [0, 0.07508265475432079, 0, 0], [0, 0, 0.07048974186182022, 0], [0.05418487141529719, 0, 0, 0], [0, 0, 0.05731044212977091, 0], [0, 0.06089730064074198, 0, 0], [0, 0.07361803700526555, 0, 0], [0, 0, 0, 0.09542109320561092], [0, 0, 0, 0.08611976976195972], [0, 0, 0.08582866564393044, 0], [0, 0, 0.06775024409095447, 0], [0.05935114622116089, 0, 0, 0], [0, 0, 0.055051843325297035, 0], [0, 0.05941355973482132, 0, 0], [0, 0, 0, 0.05964029828707377], [0, 0.051138186206420265, 0, 0], [0, 0, 0.1449235863983631, 0], [0, 0, 0.16886828218897185, 0], [0, 0, 0.17150993396838507, 0], [0, 0.07219794144233067, 0, 0], [0, 0, 0.07045306265354156, 0], [0, 0, 0, 0.08089139809211095], [0, 0.07778434952100118, 0, 0], [0, 0, 0.08693968504667282, 0], [0, 0, 0.09603794415791829, 0], [0, 0, 0.09240440651774406, 0], [0, 0, 0, 0.058710042387247086], [0, 0.059909665336211525, 0, 0], [0.058315858244895935, 0, 0, 0], [0, 0.07387932141621907, 0, 0], [0, 0.051658835262060165, 0, 0], [0, 0, 0.06555908049146335, 0], [0, 0.04428854087988535, 0, 0], [0, 0, 0.03956529373923937, 0], [0, 0, 0, 0.026596753547588985], [0, 0.033451881259679794, 0, 0], [0, 0.03880758583545685, 0, 0]]
    # for j in range(len(old_W1) - 15):
    #     for i in range(4):
    #         if old_W1[j + 15][i] < 1:
    #             old_W1[j + 15][i] = old_W1[j + 15][i] * 2

    X1 = np.array([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 1, 0, 0]])
    W1 = np.array(old_W1)




    plot_saliency(X1.T, W1.T)
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
