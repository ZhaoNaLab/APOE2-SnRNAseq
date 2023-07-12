import cv2
import numpy as np
import pandas as pd
import os
import concurrent.futures
from PIL import Image

def process_image(file):
    try:
        directory, filename = os.path.split(file)
        sample_id = os.path.split(directory)[1].split(",")[0]
        image_name = os.path.splitext(filename)[0]

        img = cv2.imread(file, cv2.IMREAD_UNCHANGED)
        if img is None:
            print(f"Error: Unable to open image file {file}")
            return None, None

        if len(img.shape) == 2:
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
        img = cv2.convertScaleAbs(img, alpha=(255.0/65535.0)) 

        image_mask = np.zeros((img.shape[0], img.shape[1], 3), dtype=np.uint8)
        num_labels, labels = cv2.connectedComponents(cv2.cvtColor(img, cv2.COLOR_BGR2GRAY))

        data = {'Sample ID': [], 'Image Name': [], 'Centroid ID': [], 'Radius': [], 'Area': [], 'Perimeter': [], 'Circularity': []}
        circle_data = {'Sample ID': [], 'Image Name': [], 'Circle ID': [], 'X': [], 'Y': [], 'Radius0': [], 'Radius1': [], 'Radius2': []}

        for label in range(1, num_labels):
            single_cluster = np.uint8(labels == label) * 255
            contours, _ = cv2.findContours(single_cluster, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            if contours:
                cnt = contours[0]
                area = cv2.contourArea(cnt)

                if area >= 11000:
                    (x, y), radius = cv2.minEnclosingCircle(cnt)
                    perimeter = cv2.arcLength(cnt, True)
                    circularity = (4 * np.pi * area) / (perimeter ** 2) if perimeter else 0

                    data['Sample ID'].append(sample_id)
                    data['Image Name'].append(image_name)
                    data['Centroid ID'].append(label)
                    data['Radius'].append(radius)
                    data['Area'].append(area)
                    data['Perimeter'].append(perimeter)
                    data['Circularity'].append(circularity)

                    circle_data['Sample ID'].append(sample_id) 
                    circle_data['Image Name'].append(image_name)
                    circle_data['Circle ID'].append(label)
                    circle_data['X'].append(int(x))
                    circle_data['Y'].append(int(y))
                    circle_data['Radius0'].append(radius)
                    circle_data['Radius1'].append((radius+40))
                    circle_data['Radius2'].append((radius+80))

                    M = cv2.moments(cnt)
                    if M["m00"] != 0:
                        cx = int(M['m10']/M['m00'])
                        cy = int(M['m01']/M['m00'])

                        center = (int(x), int(y))
                        cv2.circle(img, center, int(radius), (0, 0, 255), 2)

                        mask1 = np.zeros((img.shape[0], img.shape[1]), dtype=np.uint8)
                        mask2 = np.zeros((img.shape[0], img.shape[1]), dtype=np.uint8)
                        mask3 = np.zeros((img.shape[0], img.shape[1]), dtype=np.uint8)

                        cv2.circle(mask1, center, int(radius), 255, -1)
                        cv2.circle(mask2, center, int(radius+40), 255, -1)
                        cv2.circle(mask3, center, int(radius+80), 255, -1)

                        image_mask[:, :, 0] = cv2.bitwise_or(image_mask[:, :, 0], mask1)
                        image_mask[:, :, 1] = cv2.bitwise_or(image_mask[:, :, 1], mask2)
                        image_mask[:, :, 2] = cv2.bitwise_or(image_mask[:, :, 2], mask3)

                        cv2.circle(img, (cx, cy), 5, (0, 0, 255), -1)
                        cv2.putText(img, f"Area: {area}", (cx+10, cy+10), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (0, 0, 255), 2)
                        cv2.putText(img, f"Circularity: {circularity:.2f}", (cx+10, cy+60), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (0, 0, 255), 2)

        pil_image = Image.fromarray(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
        pil_image_mask = Image.fromarray(image_mask)

        pil_image.save(file.replace('.tiff', '_circled.tiff'), compression='tiff_deflate')
        pil_image_mask.save(file.replace('.tiff', '_mask.tiff'), compression='tiff_deflate')

        return pd.DataFrame(data), pd.DataFrame(circle_data)

    except Exception as e:
        print(f"Error processing image {file}: {str(e)}")
        return None, None

Dir = ".../IFAnalysis"

tiff_files = [os.path.join(root, file) for root, dirs, files in os.walk(Dir) for file in files if file.endswith(".tiff") and "_processedBlue" in file]

dataframes = []
mask_dataframes = []
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    results = executor.map(process_image, tiff_files)
    for result in results:
        df, circle_df = result
        if df is not None and circle_df is not None:
            dataframes.append(df)
            mask_dataframes.append(circle_df)

df_all = pd.concat(dataframes, ignore_index=True)
circle_df_all = pd.concat(mask_dataframes, ignore_index=True)

df_all.to_csv('.../Results/PlaqueROIInfo.csv', index=False)
circle_df_all.to_csv('.../Results/MaskInfo.csv', index=False)
