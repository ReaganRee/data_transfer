import os
import sys
folder = '/Users/reagan/Desktop/3UTRBERT_visualiztion/motif_stack/motifs_PWM'
for filename in os.listdir(folder):
    infilename = os.path.join(folder,filename)
    if not os.path.isfile(infilename): continue
    #oldbase = os.path.splitext(filename)
    newname = infilename.replace('.txt', '.pcm')
    output = os.rename(infilename, newname)
