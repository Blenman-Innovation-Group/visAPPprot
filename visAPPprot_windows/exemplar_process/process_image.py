from skimage.io import imread
from method1 import get_coordinates
from object_extract import get_objects
from text_extract import get_text
import os
import matplotlib.pyplot as plt
import numpy as np
import statistics

# non_pathway
# [714, 109, 740, 714, 1231, 1405, 626, 359, 785, 496, 137, 1613, 756, 86, 149, 491, 666, 1613, 947, 649, 712, 234, 821, 51]
# [235, 3, 354, 235, 996, 954, 482, 117, 364, 299, 3, 1056, 447, 3, 3, 234, 417, 1056, 664, 446, 470, 3, 519, 3]
# 689.0
# 359.0

# pathway
# [1834, 535, 1330, 560, 810, 712, 331, 1121, 837, 666, 371, 324, 518, 707, 890, 947, 763, 1859]
# [488, 269, 208, 177, 122, 273, 188, 167, 102, 270, 149, 147, 133, 355, 303, 196, 303, 341]
# 737.5
# 202.0

def main(directory):
    count = 0
    num_objects_total = []
    num_words_total = []
    img_names = []

    words_list = {} # image name : words list as comma separated string

    for filename in os.listdir(directory):
        img_path = os.path.join(directory, filename)
        # checking if it is a file

        # try:

        if os.path.isfile(img_path) and ((".png" in img_path) or (".jpg" in img_path)):                
            # shadow_img_path = '../../Images/shadows_out.png'

            # img_path = directory + filename
            print(img_path)


            shadow_img_path = img_path
            image = imread(img_path)
            shadow_img = imread(shadow_img_path)
            word_list = get_coordinates(img_path, "grayscale", 1, 1, 5)
            svg_dir = 'C:\\Users\\Kelly\\OneDrive - Yale University\\Fall2023\\CPSC490\\Foundational Code\\boundarytest\\test'
            num_objects = len(get_objects(svg_dir, image, word_list))
            num_words, words = get_text(word_list, shadow_img, img_path)
            num_objects_total.append(num_objects)
            num_words_total.append(num_words)
            # img_names.append("#" + filename.split(".")[0].split("_")[1])
            count += 1
            print(count)
            print("-----------------------")
            
            words_list[filename] = '||'.join(words) # [word, x, y]

                # print(len(word_list))

                # print(num_objects, num_words)
                # if ((num_words == 3) or (num_words >= (0.5*num_objects))):
                #     cmd = "mv " + img_path + " ../../non_pathway/" + filename
                #     os.system(cmd)
                # else:
                #     # probably more shapes than text
                #     cmd = "mv " + img_path + " ../../pathway/" + filename
                #     os.system(cmd)
                    
        # except Exception as error:
        #     print("ERROR -----------------")
        #     print(error)

    # width = 0.1
    # x = np.arange((count)) 
    # print(x)
    # plt.bar(x-0.1, num_objects_total, width, color='cyan') 
    # plt.bar(x+0.1, num_words_total, width, color='orange') 
    # plt.legend(['num_objects', 'num_words'])
    # ax = plt.axes() 
    # ax.set_xticks(list(x))
    # ax.set_xticklabels(img_names) 
    # plt.gcf().set_size_inches(10, 5)
    # plt.savefig('pdf_hist.png', dpi=200)

    # print(num_objects_total)
    # print(num_words_total)
    # print(statistics.median(num_objects_total))
    # print(statistics.median(num_words_total))


    # write words list to file
    words_list_f = open("./words_list.txt", "w")

    for filename in words_list:
        string = filename + "||" + words_list[filename] + '\n'
        words_list_f.write(string)

    words_list_f.close()

    print(words_list)

    


     

directory = "./imgs/"
main(directory)



# draw out the diagram of the full workflow
# incorporate preprocessing into this file as well (shadows removal, blurring/grayscaling)