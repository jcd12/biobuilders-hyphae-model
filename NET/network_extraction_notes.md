

python2 binarize/binarize_adaptive.py ../../../image_processing/TM_12_0_5_57_5_156.jpg -g 7 -s 1 -t 51 -m 70000 -c 30


python2 gegui/gegui.py data/results/igem_image/

python2 analyze/analyze.py data/results/igem_image/TM_12_0_5_57_5_156_graph_p3_r0.gpickle -dest data/results/igem_image/


# lige nu beregner vi vinkler for alle nabopar - det skal vi ikke 
    # vi kan fjerne naboparret, hvis vinkel er tættest på 180 grader
    # for hver af de resterende naboer, kan vi så tage vinklen til den fra naboparret, der giver den mindste vinkel, da dette most likely vil være branching angle?!? ikke helt godt 

