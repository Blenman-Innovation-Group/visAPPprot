# importing relevant libraries
import matplotlib.pyplot as plt
import numpy as np
import cv2
import pytesseract
import easyocr
from exemplar_process.util import draw_ocr_boxes, add_white_bg, grayscale, binthres, otsu, rescale
from copy import deepcopy
import torch

# necessary to prevent crashing during multiple word extraction calls
torch.set_num_threads(1)
torch.set_num_interop_threads(1)

reader = easyocr.Reader(['en'])

# global vars
# pytesseract.pytesseract.tesseract_cmd = r'C:\Users\Kelly\AppData\Local\Programs\Tesseract-OCR\tesseract.exe'

# returns bounding boxes and text predictions from EasyOCR
def easyocr_text_pred(img, hths, wths):
    coords = reader.readtext(img, canvas_size=3840, 
                             height_ths=hths, width_ths=wths)

    return coords

def check_bounds(coord, bound):
    return max(0, min(int(coord), bound))

# reshapes easyocr predictions into set format
# ensure coordinates do not exceed the height and width of image
def recons_eocr_coords(coordinates, width, height):
    reconstruct = []
    for point in coordinates:
        coords, _, _ = point
        reconstruct.append(np.array([[check_bounds(coord[0], width), 
                     check_bounds(coord[1], height)] 
                    for coord in coords], dtype=int))
    return reconstruct

def order_points(rectangle_coords):
    if isinstance(rectangle_coords, np.ndarray):
        rectangle_coords = rectangle_coords.tolist()
    # Find the three coordinates with the smallest y-axis values
    sorted_coords = sorted(rectangle_coords, key=lambda coord: coord[1])[:3]

    # Calculate distances between the three points using Pythagoras' theorem
    def distance(x1, y1, x2, y2):
        return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
    
    distances = [
        distance(sorted_coords[i][0], sorted_coords[i][1], sorted_coords[j][0], sorted_coords[j][1])
        for i in range(3)
        for j in range(i + 1, 3)
    ]
    # Find the second longest distance and return the indices of the points
    second_longest_distance = sorted(distances)[-2]
    lst = [(0,1),(1,2),(0,2)]
    indexes = (0,0)
    for x in lst:
        if distance(sorted_coords[x[0]][0], sorted_coords[x[0]][1], sorted_coords[x[1]][0], sorted_coords[x[1]][1]) == second_longest_distance:
            indexes = x
            break
        else:
            continue
        
    coordinates = sorted((sorted_coords[indexes[0]],sorted_coords[indexes[1]]), key=lambda coord: coord[0])[0]
    index = rectangle_coords.index(coordinates)
    if index == 0:
        return rectangle_coords
    else:
        # Reorder the points based on the second longest distance
        rectangle_coords = rectangle_coords[index:] + rectangle_coords[:index]
    return rectangle_coords

# performs thresholding and image inversion on the region of interest
def transform_roi(roi):
    roi_thresh = cv2.adaptiveThreshold(roi, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 21, 10)
    roi_inv = cv2.bitwise_not(roi)
    return [roi, roi_thresh, roi_inv]

def pythagorean(a, b):
    return int((a**2 + b**2)**0.5)

def warp_rois(img_list, coordinates):
    # array type must be np.float32 for the getPerspectiveTransform method
    corners = np.array(coordinates, dtype=np.float32)
    roi_width = pythagorean(corners[1][0] - corners[0][0], corners[1][1] - corners[0][1])
    roi_height = pythagorean(corners[3][0] - corners[0][0], corners[3][1] - corners[0][1])
    rectified_corners = np.array([[0, 0], [roi_width, 0], 
                                  [roi_width, roi_height], 
                                  [0, roi_height]], dtype=np.float32)
    transformation_matrix = cv2.getPerspectiveTransform(corners, rectified_corners)

    warped = []
    for img in img_list:
        warped.append(cv2.warpPerspective(img, transformation_matrix, (roi_width, roi_height)))
    return warped

def tesseract_ocr(roi_list):
    custom_config = r'--oem 1 --psm 7 -c preserve_interword_spaces=1'
    results = []
    for roi in roi_list:
        # results.append(pytesseract.image_to_data(roi, config=custom_config, output_type=pytesseract.Output.DICT))
        data = pytesseract.image_to_data(roi, config=custom_config)
        # convert string to dictionary
        data = data.strip().split("\t")

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

        for k in range(j, len(header)):
            dic[header[k]] = ""
        results.append(deepcopy(dic))

    return results

def check_angle(coordinates, img, slant, scale):
    word_list = []
    # reader = easyocr.Reader(['en'])
    for box in coordinates:
        if abs(box[0][0] - box[3][0]) < slant:
            roi = img[box[0][1]:box[2][1], box[0][0]:box[1][0]]
            if np.any(roi):
                roi_list = transform_roi(roi)
                results = tesseract_ocr(roi_list)
                conf_scores = [result['conf'] for result in results]
                index = conf_scores.index(max(conf_scores))
                border_im = add_white_bg(roi_list[0], 100)
                coords = reader.readtext(border_im, height_ths=5, width_ths=5)

                if coords and (not all(v == "" for v in conf_scores)):
                    if float(max(conf_scores)) > float(coords[0][2]) * 100:
                        word_list.append([(results[index]['text']).strip(), box/scale, "not slanted"])
                    else:
                        word_list.append([coords[0][1].strip(), box/scale, "not slanted"])

        else:
            points = order_points(box)

            roi_list = transform_roi(img)
            

            warped_rois = warp_rois(roi_list, points)
            results = tesseract_ocr(warped_rois)

            conf_scores_all = [result['conf'] for result in results]
            conf_scores = []
            for score in conf_scores_all:
                if len(score) > 0:
                    conf_scores.append(score)

            if len(conf_scores) == 0:
                continue

            index = conf_scores.index(max(conf_scores)) # finding the max confidence score of tesseract results
            border_im = add_white_bg(warped_rois[0], 100) # extending the border of the ROI
            cv2.imwrite("roi-extended.jpg", border_im)
            coords = reader.readtext(border_im) # EasyOCR applied on just the ROI

            if coords: # EasyOCR often detects the lack of text better than Tesseract
                if float(max(conf_scores)) > (float(coords[0][2]) * 100) : 
                    word_list.append([(results[index]['text']).strip(), np.divide(points, scale), "slanted"])
                else:
                    word_list.append([coords[0][1].strip(), np.divide(points, scale), "slanted"])
    return word_list

def run_filter(img, filt):
    if filt == "grayscale":
        filt_img = grayscale(img)
    elif filt == "binthres":
        filt_img = binthres(img)
    elif filt == "otsu":
        filt_img = otsu(img)
    else:
        filt_img = img
    return filt_img

# def get_scale(image):
#     height, width = image.shape

def get_coordinates(img_path, filt, hths, wths, slant):
    scale = 1

    im = cv2.imread(img_path)
    im = rescale(im, scale, scale)
    filt_img = run_filter(im, filt)
    height, width, _ = im.shape

    eocr_pred_set = easyocr_text_pred(filt_img, hths, wths)
    coords_eocr = recons_eocr_coords(eocr_pred_set, width, height)
    draw_ocr_boxes(coords_eocr, im)

    final_coords = check_angle(coords_eocr, filt_img, slant, scale)    

    return final_coords

# word_list = get_coordinates('./Images/354434.jpg', "grayscale", 1, 1, 5)