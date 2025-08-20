import cv2
import numpy as np
from skimage.io import imread
import matplotlib.pyplot as plt
import csv
import pytesseract
import easyocr
import os
from copy import deepcopy
from exemplar_process.method1 import reader

# https://stackoverflow.com/questions/58613825/get-boundary-from-canny-edges-and-remove-the-background-of-an-image
# https://stackoverflow.com/questions/43108751/convert-contour-paths-to-svg-paths
# pytesseract.pytesseract.tesseract_cmd = r'C:\Users\Kelly\AppData\Local\Programs\Tesseract-OCR\tesseract.exe'
custom_config = r'--oem 1 --psm 7 -c preserve_interword_spaces=1 -c load_system_dawg=0 -c load_freq_dawg=0 -c load_punc_dawg=0 -c load_number_dawg=0 -c language_model_penalty_non_freq_dict_word=0 -c language_model_penalty_non_dict_word=0 -c language_model_penalty_punc=0 language_model_penalty_case=0 -c tessedit_char_blacklist="|" -l ell+eng '
# easyocr_reader = easyocr.Reader(['en'])


def easyocr_box_pred(img, hths, wths):
    coords = easyocr_reader.detect(img, canvas_size=3840, 
                             height_ths=hths, width_ths=wths)
    return coords

empty_test = {}
def check_bounds(wordbounds, x, y, w, h):
    for i, box in enumerate(wordbounds):
        topleft, _, botright, _ = box[1]
        if (int(topleft[0]) < x and int(topleft[1]) < y) and (int(botright[0]) > x + w and int(botright[1]) > y + h):
            if i in empty_test:
                coords = list(empty_test[i])
                if x < coords[0]:
                    coords[0] = x
                if y < coords[1]:
                    coords[1] = y
                if x+w > coords[2]:
                    coords[2] = x+w
                if y+h > coords[3]:
                    coords[3] = y+h
                empty_test[i] = tuple(coords)
            else:
                empty_test[i] = (x, y, x+w, y+h)
            return 'letter'
    return 'object'

def get_text(word_list, image, og_path):
    global empty_test # global vs local definitions think about it
    empty_test = {}

    og_img = cv2.cvtColor(imread(og_path), cv2.COLOR_BGR2RGB)

    img = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # applying otsu threshold to the image
    blur = cv2.bilateralFilter(gray, 3, 10, 10)
    thresh_gray = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    thresh_blur = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    
    # cv2.imshow("blur", blur)
    # cv2.imshow("thresgray", thresh_gray)
    # cv2.imshow("thresblur", thresh_blur)
    # cv2.waitKey(0)

    # using openCV to find contours 
    cnts, hierarchy = cv2.findContours(thresh_blur, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    hierarchy = hierarchy[0]
    zip_contours = zip(cnts, hierarchy)
    
    for i, node in enumerate(zip_contours, 0):
        c = node[0]

        area = cv2.contourArea(c)
        rx, ry, rw, rh = cv2.boundingRect(c)
        cnt_type = check_bounds(word_list, rx, ry, rw, rh)
        # print(rx, ry, rw, rh)
        # cv2.rectangle(og_img,(rx,ry),(rx+rw,ry+rh),(255,0,0),2)
        # cv2.imshow('rect', og_img)
        # cv2.waitKey(0)
        if cnt_type == "letter":
            mask = np.zeros(img.shape, dtype=np.uint8)
            cv2.drawContours(mask, [c], -1, (255,255,255))
            result = np.bitwise_and(mask, img)
            # cv2.imshow('letter', result[ry-500:ry+rh+500,rx-500:rx+rw+500])
            # cv2.imshow('letter strict', result[ry:ry+rh,rx:rx+rw])
            # cv2.waitKey(0)
    
    word_canvas = np.zeros(img.shape, dtype=np.uint8)
    word_canvas.fill(255)
    
    x_buff = (img.shape[0] // 720) * 2
    y_buff = (img.shape[1] // 720) * 2
    text_size = (img.shape[0] // 720 * 0.5)

    words = [] # [word, x, y]

    for i, item in enumerate(empty_test.values()):
        roi = og_img[item[1]-y_buff:item[3]+y_buff, item[0]-x_buff:item[2]+x_buff]
        
        # result = pytesseract.image_to_data(roi, config=custom_config, output_type=pytesseract.Output.DICT)
        result = pytesseract.image_to_data(roi, config=custom_config)

        data = result.strip().split("\t")
        # print(data)

        header = ["level", "page_num", "block_num", "par_num", "line_num", "word_num", "left", "top", "width",   "height", "conf", "text"]
        results = []
        dic = {}
        j = 0
        for i in range(len(header), len(data)):
            val = data[i]
            if "\n" in val:
                # newline in return data
                val = val.split("\n")
                val = list(filter(None, val))

                # fill the rest of the dict with empty 
                if len(val) > 1:
                    dic[header[j]] = val[0]
                    for k in range(j+1, len(header)):
                        dic[header[k]] = ""
                else:
                    for k in range(j, len(header)):
                        dic[header[k]] = ""

                results.append(deepcopy(dic))
                j = 0
                if len(val) > 1:
                    dic[header[j]] = val[1]
                else:
                    dic[header[j]] = val[0]
                j += 1
            else:
                dic[header[j]] = val
                j += 1

        results.append(deepcopy(dic))

        word = ""

        for result in results:
            if "text" in result:
                word = word + " " + result['text']

        if len(word.strip()) > 0:
            word_x = (item[0] + item[2])/2
            word_y = (item[1] + item[3])/2
            words.append(word.strip())
            print(word.strip())
            words.append(str(word_x))
            words.append(str(word_y))

        # cv2.imshow('roi' + str(i), roi)
        # cv2.waitKey(0)
        # cv2.putText(word_canvas, ' '.join(result['text']).strip(), (item[0], item[3]), cv2.FONT_HERSHEY_SIMPLEX, text_size, (0, 0, 0), 1)

    # cv2.imshow('words', word_canvas)
    # # cv2.imwrite('C:\\Users\\Kelly\\OneDrive - Yale University\\Fall2023\\CPSC490\\Foundational Code\\Images\\wordcanvastest.jpg', word_canvas)
    # cv2.waitKey(0)

    # get number of letters in total
    num_letters = 0
    for i in range(0, len(words), 3):
    # for line in words:
        words_line = words[i]
        words_line = words_line.strip().split(" ")
        for word in words_line:
            num_letters += len(word)

    return num_letters, words

# get_text()  

# TODO list
# get center coordinate of text
    # how should text be presented? within spreadsheet with coordinates?
# remove shadowing from images... (incorporate into workflow)
    # https://stackoverflow.com/questions/44752240/how-to-remove-shadow-from-scanned-images-using-opencv
# (ADJUST THIS) set up bounds based on the shape (and resolution?) of the image
# TODO NEXT STEP: slanted text is a problem
# output central coordinates of text 
# meeting goal: have 5 separate images that can be run through one method and each get solid outputs in one go

'''
1. read original image and image processed for shadows

shadow processed image is read as input, rois are obtained from the original image


'''

# how to get the og colors? do i need the og colors? or save it as a separate file?

# run through program and see how workflow goes
# watch the video 