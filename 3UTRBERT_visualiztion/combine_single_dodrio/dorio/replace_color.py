from PIL import Image
import numpy as np

img = Image.open(
    '/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/bertviz_single/101_6_2.png')
arr = np.asarray(img)
# 打印数组的形状
print('数组的维度尺寸为：')
print(arr.shape)
print('\n')
print('打印像素值')
for i in arr:
    for j in i:
        if j[0] == 0 and j[1] == 0 and j[2] == 0:
            j[0] = 255
            j[1] = 255
            j[2] = 255
#plt.imsave('E:/1235.png',arr1)
data = Image.fromarray(arr)
data.save('/Users/reagan/Desktop/3UTRBERT_visualiztion/combine_single_dodrio/bertviz_single/101_6_2_replaced.png')
