import os, glob

filename = "out"
os.chdir("exrFiles")

for file in glob.glob("*.exr"):
	print(file[0:len(file)-4])
	os.system(f"./../magick {file} -auto-gamma ../{file[0:len(file)-4]}.png")


os.chdir("../")
os.system(f"ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4")