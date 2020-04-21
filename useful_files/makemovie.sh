#./ffmpeg -r 2 -i %1d.png -pix_fmt yuv420p -r 10 dfrac.mp4
./ffmpeg -r 2 -i %1d.png -f mp4 -q:v 0 -vcodec mpeg4 -vb 20M -r 10 dfrac.mp4
