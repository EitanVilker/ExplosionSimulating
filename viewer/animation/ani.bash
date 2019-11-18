#!/bin/bash
$ chmod +x ani.bash
echo "Create SLAX animation!"
echo "Param 1: $1"
ffmpeg -framerate 50 -i $1/_images/%4d.png -c:v libx264 -vf scale=640:-1 -pix_fmt yuv420p $1/_images/_video.mp4
