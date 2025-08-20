import cv2
import numpy as np
from skimage.io import imread
import matplotlib.pyplot as plt
import csv
import os

# https://stackoverflow.com/questions/58613825/get-boundary-from-canny-edges-and-remove-the-background-of-an-image
# https://stackoverflow.com/questions/43108751/convert-contour-paths-to-svg-paths

def write_csv(write_data):
    with open('C:\\Users\\Kelly\\OneDrive - Yale University\\Fall2023\\CPSC490\\Foundational Code\\boundarytest\\test.csv', 'w', newline='') as f_w:
        headers = ['filename', 'x-coord', 'y-coord', 'width, height', 'corners (clockwise): tl, tr, br, bl']
        writer = csv.DictWriter(f_w, fieldnames=headers)
        writer.writeheader()
        writer.writerows(write_data)

def check_bounds(wordbounds, x, y, w, h):
    for box in wordbounds:
        topleft, _, botright, _ = box[1]
        if (int(topleft[0]) < x and int(topleft[1]) < y) and (int(botright[0]) > x + w and int(botright[1]) > y + h):
            return 'letter'
    return 'object'

def write_svg(contour, off_x, off_y):
    # rx, ry, rw, rh = cv2.boundingRect(contour)
    svgpath = '<path d="M'
    for j in range(len(contour)):
        x, y = contour[j][0]
        svgpath += f"{x-off_x} {y-off_y} "
    svgpath += '" style="stroke:black"/>'
    return svgpath

def dfs_family(word_list, visited, cnt_list, h_list, nodeidx, svgpaths, off_x, off_y):
    if nodeidx not in visited:
        # toplevel node
        print(nodeidx)
        visited.add(nodeidx)
        contour = cnt_list[nodeidx]
        hierarchy = h_list[nodeidx]
        svgpaths.append(write_svg(contour, off_x, off_y))
        # child nodes
        child = hierarchy[2]
        # TODO: reiterate over code, make less redundant
        while(child != -1):
            child_contour = cnt_list[child]
            area = cv2.contourArea(child_contour)
            rx, ry, rw, rh = cv2.boundingRect(child_contour)
            cnt_type = check_bounds(word_list, rx, ry, rw, rh)
            if area > 0 and cnt_type != 'letter':
                # svgpaths.append(write_svg(child_contour, off_x, off_y))
                # search for grandchildren
                grandchild = h_list[child][2]
                if grandchild != -1:
                    svgpaths = dfs_family(word_list, visited, cnt_list, h_list, grandchild, svgpaths, off_x, off_y)
            # otherwise go to next child of the parent
            child = h_list[child][0]
    return svgpaths

def get_objects(svg_writedir, image, word_list):

    img = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    why = np.zeros(img.shape, dtype=np.uint8)
    # applying otsu threshold to the image
    blur = cv2.bilateralFilter(gray, 3, 10, 10)
    # thresh_gray = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    thresh_blur = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    
    # cv2.imshow("blur", blur)
    # # cv2.imshow("thresgray", thresh_gray)
    # cv2.imshow("thresblur", thresh_blur)
    # cv2.waitKey(0)


    contour_info = []
    # using openCV to find contours 
    cnts, hierarchy = cv2.findContours(thresh_blur, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    hierarchy = hierarchy[0]
    zip_contours = zip(cnts, hierarchy)
    for i, node in enumerate(zip_contours, 0):
        # creating new svg file for each contour found within the image
        c = node[0]
        
        # use hierarchy to combine objects together into singular svg
        h = node[1]
        write_contour = {}

        area = cv2.contourArea(c)
        rx, ry, rw, rh = cv2.boundingRect(c)
        cnt_type = check_bounds(word_list, rx, ry, rw, rh)

        if h[3] == -1 and cnt_type == 'object' and area > 0:
            write_contour = {}
            write_contour['width, height'] = (rw, rh)
            write_contour['corners (clockwise): tl, tr, br, bl'] = [(rx, ry), (rx+rw, ry), (rx+rw, ry+rh), (rx, ry+rh)]
            filename = os.path.join(svg_writedir, 'object' + str(i) + '.svg')
            write_contour['filename'] = filename

            # with open(filename, "w+") as f:

            #     # svg file header
            #     f.write(f'<svg fill="none" width="{rw}" height="{rh}">')
            #     f.write('<path d="M')
            #     for j in range(len(c)):
            #         x, y = c[j][0]
            #         f.write(f"{x-rx} {y-ry} ")
            #     f.write('" style="stroke:black"/>')
            #     visited = set()
            #     svgList = []
            #     svg_paths = dfs_family(word_list, visited, cnts, hierarchy, i, svgList, rx, ry)
                
            #     for svgcontour in svg_paths[1:]:
            #         f.write(svgcontour)
                
            #     f.write("</svg>")
                # equation for finding center coordinates of contour (https://www.geeksforgeeks.org/python-opencv-find-center-of-contour/)
                # find center coordinates for each contour and adds them to csv file 
            moment = cv2.moments(c)
            if moment['m00'] != 0:
                cx = int(moment['m10']/moment['m00'])
                cy = int(moment['m01']/moment['m00'])
                write_contour['x-coord'] = cx
                write_contour['y-coord'] = cy
            else:
                write_contour['x-coord'] = 'n/a'
                write_contour['y-coord'] = 'n/a'
            # for visual confirmation, creates mask based off the contour 
            mask = np.zeros(img.shape, dtype=np.uint8)
            cv2.drawContours(mask, [c], -1, (255,255,255))
            # cv2.drawContours(img, [c], -1, (255,255,255), -1)

            result = np.bitwise_and(mask, img)
            # cv2.imshow('mask', result)
            # cv2.waitKey(0)

            contour_info.append(write_contour)
            cv2.drawContours(why, [c], -1, (255,255,255))
    # write_csv(contour_info)

    # cv2.imshow('final', why)
    # invert = cv2.bitwise_not(why)
    # cv2.imshow('inversion', invert)
    # # cv2.imwrite('C:\\Users\\Kelly\\OneDrive - Yale University\\Fall2023\\CPSC490\\Foundational Code\\Images\\outputcontours-test.jpg', invert)
    # cv2.waitKey(0)

    return cnts
# get('C:\\Users\\Kelly\\OneDrive - Yale University\\Fall2023\\CPSC490\\Foundational Code\\boundarytest\\354434')  

# TODO list
# get center coordinate of text
    # how should text be presented? within spreadsheet with coordinates?
# remove shadowing from images...
    # https://stackoverflow.com/questions/44752240/how-to-remove-shadow-from-scanned-images-using-opencv
# set up a cleaner workflow for all the steps

# find header that's missing to view svg directly
# transfer files to yale drive


# meeting goal: have 5 separate images that can be run through one method and each get solid outputs in one go
'''
1. read original image and image processed for shadows

shadow processed image is read as input, rois are obtained from the original image


'''

# how to get the og colors? do i need the og colors? or save it as a separate file?

# run through program and see how workflow goes
# watch the video 