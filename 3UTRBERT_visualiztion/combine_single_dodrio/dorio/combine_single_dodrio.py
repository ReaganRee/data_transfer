# Python实现将一张图片放到另一张图片指定的位置上并合成一张图
import random

from PIL import Image
import numpy as np
random.seed(1)
Image.MAX_IMAGE_PIXELS = None
background = Image.new('RGB', (30000, 30000), (255, 255, 255))
background.save(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/background.png")
importance_score = np.loadtxt(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/importance_score_seq1.txt")
print(importance_score)

background_empty = Image.open(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/background.png")
for i in range(12):
    for j in range(12):
        img_to_paste = Image.open(
            "/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/Layer_{layer}/head_{head}/crop.jpg".format(
                layer=i, head=j))
        if importance_score[i][j] > 0.99:
            factor = random.uniform(0.8, 1)
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * factor), int(2500 * importance_score[i][j] * factor)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * factor),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * factor))

        if importance_score[i][j] >= 0.5:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j]), int(2500 * importance_score[i][j])))
            position = (int(1250 + j * 2500 - 1250 * importance_score[i][j]),
                        int(1250 + i * 2500 - 1250 * importance_score[i][j]))
        if importance_score[i][j] < 0.5:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 1.5), int(2500 * importance_score[i][j] * 1.5)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 1.5),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 1.5))

        if importance_score[i][j] < 0.3:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 2.2), int(2500 * importance_score[i][j] * 2.2)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 2.2),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 2.2))

        if importance_score[i][j] < 0.2:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 3.5), int(2500 * importance_score[i][j] * 3.5)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 3.5),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 3.5))


        if importance_score[i][j] < 0.15:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 5), int(2500 * importance_score[i][j] * 5)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 5),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 5))
        if importance_score[i][j] < 0.1:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 7), int(2500 * importance_score[i][j] * 7)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 7),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 7))

        if importance_score[i][j] < 0.05:
            img_to_paste = img_to_paste.resize((int(2500 * importance_score[i][
                j] * 9), int(2500 * importance_score[i][j] * 9)))
            position = (
                int(1250 + j * 2500 - 1250 * importance_score[i][j] * 9),
                int(1250 + i * 2500 - 1250 * importance_score[i][j] * 9))


        background_empty.paste(img_to_paste, position, mask=None)

background_empty.save(
    "/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/Dodrio_completed.png")  # 合成后的图片路径以及文件名

# 前两个坐标点是左上角坐标
# 后两个坐标点是右下角坐标
# width在前， height在后
# box = (625, 0, 3125, 2500)
# for i in range(12):
#     for j in range(12):
#         img_crop = Image.open("/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq2/Layer_{layer}/head_{head}/seq2_{layer}_{head}.png".format(layer=i, head=j))
#         img_crop = img_crop.convert("RGB")
#         region = img_crop.crop(box)
#         region.save("/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq2/Layer_{layer}/head_{head}/crop.jpg".format(layer=i, head=j))

# img_crop = Image.open("/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/Layer_10/head_5/seq1_10_5.png")
# print(img_crop.size)
# img_crop = img_crop.resize((3750,2500)) # 转化图片
# img_crop = img_crop.convert("RGB")
# region = img_crop.crop(box)
# region.save("/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/dorio/seq1/Layer_10/head_5/crop.jpg")
