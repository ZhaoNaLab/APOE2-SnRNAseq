import cv2
import numpy as np
import pandas as pd
import os
import re

def generate_mask(image_shape, x, y, radius):
    mask = np.zeros(image_shape, dtype=np.uint8)
    cv2.circle(mask, (x, y), int(radius), 255, -1)
    return mask

def generate_outside_mask(image_shape, mask):
    return cv2.bitwise_not(mask)

Dir = '.../IFAnalysis'  

tiff_files = [os.path.join(root, file) for root, dirs, files in os.walk(Dir) for file in files if re.search(r'\d+\.tiff$', file)]
Sample_ID = [os.path.basename(os.path.dirname(file)).split(",")[0] for file in tiff_files]
Image_Name = [os.path.splitext(os.path.basename(file))[0] for file in tiff_files]

maskInfo = pd.read_csv('.../Results/MaskInfo.csv')
maskInfo = pd.melt(maskInfo, id_vars=['Sample ID', 'Image Name', 'Circle ID', 'X', 'Y'], value_vars=['Radius0', 'Radius1', 'Radius2'], var_name='RadiusGroup', value_name='Radius')
maskInfo["Image Name"] = maskInfo["Image Name"].str.replace("_processedBlue", "")

results_inside= []
results_outside=[]

for sample_id in maskInfo['Sample ID'].unique():
    sample_maskInfo = maskInfo[maskInfo['Sample ID'] == sample_id]

    for imageName, imageGroup in sample_maskInfo.groupby('Image Name'):
        tiff_file = next((tf for tf in tiff_files if os.path.splitext(os.path.basename(tf))[0] == imageName), None)
        if tiff_file is None:
            continue

        img = cv2.imread(tiff_file, cv2.IMREAD_UNCHANGED)
        if img is None:
            print(f"Error reading file {tiff_file}")
            continue

        img_blue = img.copy()[:,:,0]
        img_green = img.copy()[:,:,1]
        img_red = img.copy()[:,:,2]

        thres_blue = 18849
        thres_green = 17497
        thres_red = 10486

        img_blue[img_blue < thres_blue] = 0
        img_green[img_green < thres_green] = 0
        img_red[img_red < thres_red] = 0

        image_shape = img_green.shape

        for radiusGroup, radiusGroupMaskInfo in imageGroup.groupby('RadiusGroup'):
            merged_mask = np.zeros(image_shape, dtype=np.uint8)

            for _, row in radiusGroupMaskInfo.iterrows():
                circleID, x, y, radius = row['Circle ID'], row['X'], row['Y'], row['Radius']

                if pd.isna(x) or pd.isna(y) or pd.isna(radius):
                    print(f"Skipping due to NaN values: x={x}, y={y}, radius={radius}")
                    continue
                x, y = map(int, (x, y))

                circle_mask = generate_mask(image_shape, x, y, radius)
                merged_mask = cv2.bitwise_or(merged_mask, circle_mask)

                img_blue_masked = cv2.bitwise_and(img_blue, img_blue, mask=circle_mask)
                img_green_masked = cv2.bitwise_and(img_green, img_green, mask=circle_mask)
                img_red_masked = cv2.bitwise_and(img_red, img_red, mask=circle_mask)

                green_red_overlap_mask = cv2.bitwise_and(circle_mask, ((img_red > 0) & (img_green > 0)).astype(np.uint8))
                img_green_red_masked = cv2.bitwise_and(img_green, img_green, mask=green_red_overlap_mask)

                total_blue = np.sum(img_blue_masked)
                total_green = np.sum(img_green_masked)
                total_red = np.sum(img_red_masked)
                total_green_red_overlap = np.sum(img_green_red_masked)

                count_blue = np.count_nonzero(img_blue_masked)
                count_green = np.count_nonzero(img_green_masked)
                count_red = np.count_nonzero(img_red_masked)
                count_green_red_overlap = np.count_nonzero(img_green_red_masked)

                results_inside.append([sample_id, imageName, circleID, radiusGroup, total_blue, count_blue, total_green, count_green, total_red, count_red, total_green_red_overlap, count_green_red_overlap, np.pi * radius**2])

            outside_mask = generate_outside_mask(image_shape, merged_mask)
            outside_area = np.count_nonzero(outside_mask)

            img_blue_masked_outside = cv2.bitwise_and(img_blue, img_blue, mask=outside_mask)
            img_green_masked_outside = cv2.bitwise_and(img_green, img_green, mask=outside_mask)
            img_red_masked_outside = cv2.bitwise_and(img_red, img_red, mask=outside_mask)

            green_red_overlap_mask_outside = cv2.bitwise_and(outside_mask, ((img_red > 0) & (img_green > 0)).astype(np.uint8))
            img_green_red_masked_outside = cv2.bitwise_and(img_green, img_green, mask=green_red_overlap_mask_outside)

            total_blue_outside = np.sum(img_blue_masked_outside)
            total_green_outside = np.sum(img_green_masked_outside)
            total_red_outside = np.sum(img_red_masked_outside)
            total_green_red_overlap_outside = np.sum(img_green_red_masked_outside)

            count_blue_outside = np.count_nonzero(img_blue_masked_outside)
            count_green_outside = np.count_nonzero(img_green_masked_outside)
            count_red_outside = np.count_nonzero(img_red_masked_outside)
            count_green_red_overlap_outside = np.count_nonzero(img_green_red_masked_outside)

            results_outside.append([sample_id, imageName, radiusGroup, total_blue_outside, count_blue_outside, total_green_outside, 
                                 count_green_outside, total_red_outside, count_red_outside, total_green_red_overlap_outside, count_green_red_overlap_outside, outside_area])

results_inside_df = pd.DataFrame(results_inside, columns=['SampleID', 'ImageName', 'CircleID', 'RadiusGroup', 'BlueIntensityROI', 'BlueCountROI', 'GreenIntensityROI', 
                                                      'GreenCountROI', 'RedIntensityROI', 'RedCountROI', 'green_redIntensityROI', 'green_redCountROI', 'ROIArea'])

results_outside_df = pd.DataFrame(results_outside, columns=['SampleID', 'ImageName', 'RadiusGroup', 'BlueIntensityBackground', 'BlueCountBackground', 'GreenIntensityBackground', 'GreenCountBackground', 
                                                        'RedIntensityBackground', 'RedCountBackground', 'green_redIntensityBackground', 'green_redCountBackground', "BackgroundArea"])

## Merge the ROI and background data
df_merge = results_inside_df.merge(results_outside_df, on=['SampleID', 'ImageName', 'RadiusGroup'], how='outer')
df_merge.to_csv('.../Results/PlaqueAssociatedMicroglia.csv', index=False)
