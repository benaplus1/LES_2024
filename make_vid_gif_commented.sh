#!/bin/bash

rm -f output.gif output_slow.gif video.mp4

#CHANGE THIS to the folder and prefix of your output images

FOLDER_NAME=/moonbow/ascheb/les/PresFigs/3dplots_les/nolake/div3d_hires_nolake


#Create the color palette for the animated GIF

ffmpeg -y -pattern_type glob -i $FOLDER_NAME"*.png" -vf scale=1080:-1:flags=lanczos,palettegen palette.png

#Create the animated gif.

ffmpeg -pattern_type glob -i $FOLDER_NAME"*.png" -i palette.png -t 60 -filter_complex "scale=1920:-1:flags=lanczos[x];[x][1:v]paletteuse" output.gif

#CHANGE THIS if you want the gif to animate slower or faster.

convert -delay 40x100 output.gif video.gif

echo "mv video.gif $FOLDER_NAME""_video.gif"

mv video.gif "$FOLDER_NAME""_video.gif"



#Create the video. CHANGE THIS if you want to change the bitrate/frame rate of your output video.

ffmpeg -framerate 6 -pattern_type glob -i $FOLDER_NAME"*.png" -pix_fmt yuv420p -b 2500k -vcodec libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" video.mp4



mv video.mp4 "$FOLDER_NAME""_video.mp4"



echo "mv video.mp4 ""$FOLDER_NAME""_video.mp4"