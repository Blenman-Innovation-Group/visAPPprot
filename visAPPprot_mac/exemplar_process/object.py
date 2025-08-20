import cv2
import numpy as np
from skimage.io import imshow, imread
from skimage.color import rgb2gray
import matplotlib.pyplot as plt
import diplib as dip

# # canny edge detection
# thresh = cv2.Canny(gray,250,200)
# cv2.imshow("thresh", thresh)
# cv2.waitKey(0)

# sobel edge detection
# src = cv2.GaussianBlur(image, (3, 3), 0)
# gray = cv2.cvtColor(src, cv2.COLOR_BGR2GRAY)
# scale = 1
# delta = 0
# ddepth = cv2.CV_16S    
# grad_x = cv2.Sobel(gray, ddepth, 1, 0, ksize=3, scale=scale, delta=delta, borderType=cv2.BORDER_DEFAULT)
# # Gradient-Y
# # grad_y = cv2.Scharr(gray,ddepth,0,1)
# grad_y = cv2.Sobel(gray, ddepth, 0, 1, ksize=3, scale=scale, delta=delta, borderType=cv2.BORDER_DEFAULT)
# abs_grad_x = cv2.convertScaleAbs(grad_x)
# abs_grad_y = cv2.convertScaleAbs(grad_y)
# thresh = cv2.addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0)
# cv2.imshow('thresh', thresh)
# cv2.waitKey()

def grayscale(img_path):
    sample = imread(img_path)
    sample_g = rgb2gray(sample)
    fig, ax = plt.subplots(1,2,figsize=(10,5))
    ax[0].imshow(sample)
    ax[1].imshow(sample_g,cmap='gray')
    ax[0].set_title('Colored Image',fontsize=15)
    ax[1].set_title('Grayscale Image',fontsize=15)
    plt.show()
    return sample_g

def grayscale_cv(img_path):
    image = cv2.imread(img_path) 
    print(image.shape)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY) 
    return gray

def binarize(sample_g):
    fig, ax = plt.subplots(1,3,figsize=(15,5))
    sample_b = sample_g > 0.6
    ax[0].set_title('Grayscale Image',fontsize=20)
    ax[0].imshow(sample_g,cmap='gray')
    ax[1].plot(sample_g[600])
    ax[1].set_ylabel('Pixel Value')
    ax[1].set_xlabel('Width of Picture')
    ax[1].set_title('Plot of 1 Line',fontsize=15)
    ax[2].set_title('Binarized Image',fontsize=15)
    ax[2].imshow(sample_b,cmap='gray')
    plt.show()

def sobel_edges(img_path):
    image = cv2.imread(img_path) 
    src = cv2.GaussianBlur(image, (3, 3), 0)
    gray = cv2.cvtColor(src, cv2.COLOR_BGR2GRAY)

    window_name = ('Sobel Demo - Simple Edge Detector')
    scale = 1
    delta = 0
    ddepth = cv2.CV_16S    
    
    grad_x = cv2.Sobel(gray, ddepth, 1, 0, ksize=3, scale=scale, delta=delta, borderType=cv2.BORDER_DEFAULT)
    # Gradient-Y
    # grad_y = cv2.Scharr(gray,ddepth,0,1)
    grad_y = cv2.Sobel(gray, ddepth, 0, 1, ksize=3, scale=scale, delta=delta, borderType=cv2.BORDER_DEFAULT)
    
    
    abs_grad_x = cv2.convertScaleAbs(grad_x)
    abs_grad_y = cv2.convertScaleAbs(grad_y)
    
    
    grad = cv2.addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0)
    
    # cv2.imwrite('./output/sobel-edges.jpg', grad)
    # cv2.imshow(window_name, grad)


def canny_edges(img_path):
    image = cv2.imread(img_path) 
    img = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    canny = cv2.Canny(img, 30, 120)
    # cv2.imshow('canny', canny)
    # cv2.waitKey(0)
    # cv2.imwrite('./output/canny-edges.jpg', canny)

def tophat_edges(img_path):
    image = cv2.imread(img_path) 
    img = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    a = dip.Image(img)
    # a.Show()
    lines = dip.Tophat(a, 30, polarity='black')

    dip.SetBorder(lines, [0], [2])

    lines = dip.PathOpening(lines, length=100, polarity='opening', mode={'robust'})
    lines = dip.Threshold(lines, method='otsu')[0]
    dip.viewer.ShowModal(lines)
    # dip.ImageWriteJPEG(lines, './output/tophatfilter-edges.jpg')
    # print(lines)
    # dip.Show("test", lines)
    # cv2.waitKey(0)

    # cv2.imshow("Edge map", edged) 
    # cv2.waitKey(0)

    # kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(9,9))
    # dilated = cv2.dilate(img, kernel)
    # eroded=cv2.erode(dilated,kernel)
    # cv2.imshow("eroded", eroded)
    # cv2.waitKey(0)
    # cnts, _ = cv2.findContours(eroded.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    # print(cnts)
    # for c in cnts:
    #     peri = cv2.arcLength(c, True)
    #     approx = cv2.approxPolyDP(c, 0.3 * peri, True)
    #     cv2.drawContours(img, [cnts[0]], -1, (0, 255, 0), 2)
    # cv2.imshow("test", img)
    # cv2.waitKey(0)
  
    # Finding Contours 
    # Use a copy of the image e.g. edged.copy() 
    # since findContours alters the image 
    # contours, hierarchy = cv2.findContours(edged,  
    #     cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE) 
    # cv2.imshow('Canny Edges After Contouring', edged) 
    # cv2.waitKey(0) 
  
    # print("Number of Contours found = " + str(len(contours))) 
    
    # Draw all contours 
    # -1 signifies drawing all contours 
    # im_test = np.zeros((6666, 3364, 1), np.uint8)
    # cv2.drawContours(image, contours, -1, (0, 255, 0), 3) 
    
    # imS = cv2.resize(image, (960, 540)) 
    # cv2.imshow('Contours', imS) 
    # cv2.waitKey(0) 
    # cv2.destroyAllWindows() 
    
    # im2 = np.zeros((10000, 10000, 1), np.uint8)
    # for contour in contours:
    #     cv2.drawContours(im2,[contour], 0, (255, 255, 255),3)
    
    # imS = cv2.resize(im2, (1000, 1000)) 
    # cv2.imwrite("test_edge.jpg", imS)
    # cv2.imshow('Test Contour', imS) 
    # cv2.waitKey(0) 
    # print(contour)


def main():
    img_path = '../Images/354434.jpg'
    gs_img = grayscale_cv(img_path)
    # binarize(gs_img)
    sobel_edges(img_path)
    canny_edges(img_path)
    # tophat_edges(img_path)
    return

if __name__ == "__main__":
    main()