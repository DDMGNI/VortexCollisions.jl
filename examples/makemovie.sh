#!/bin/bash

img="FokkerPlanck_sinx4sinyetc_dt1E-3_nt10000_movie_"
out="FokkerPlanck_sinx4sinyetc_dt1E-3_nt10000.mp4"

opt="vbitrate=5000000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq"

mencoder -ovc lavc -lavcopts threads=4:vcodec=mpeg4:vpass=1:$opt -lavdopts threads=4 -mf w=1920:h=1200:fps=60:type=png -nosound -o /dev/null  mf://$img*.png
mencoder -ovc lavc -lavcopts threads=4:vcodec=mpeg4:vpass=2:$opt -lavdopts threads=4 -mf w=1920:h=1200:fps=60:type=png -nosound -o $out       mf://$img*.png
