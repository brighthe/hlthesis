import os
from PIL import Image

# 创建一个包含所有要合并的PNG图像文件名的列表
png_file_list = [s[:-1] for s in os.popen('ls *.png')]

# 设置GIF文件的输出文件名
output_gif = 'output.gif'

# 打开PNG图像并将它们添加到一个列表中
images = []
for png_file in png_file_list:
    img = Image.open(png_file)
    images.append(img)

# 将图像列表保存为GIF
images[0].save(output_gif, save_all=True, append_images=images[1:],
        duration=100, loop=0, transparency=0)

print(f'GIF文件 {output_gif} 已创建成功!')
