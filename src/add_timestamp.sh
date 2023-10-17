#!/bin/bash

# Script to add timestamp and frame numbers to video files

WD=${1:?Provide a working directory to process}

cd "$WD"

ls *.mp4  | cut -f 1 -d "." > temp.out

numlines=$(grep -c "^" temp.out)

for i in `seq 1 $numlines`; do
  a=$(awk -v "var=$i" 'NR==var' temp.out);
  ffmpeg -i "$a".mp4 -vf "drawtext=text='Frame %{n}, Time %{pts\:hms}': x=(w-tw)/2: y=h-(2*lh): fontcolor=white: fontsize=48" -y "$a"l.mp4
done

rm temp.out