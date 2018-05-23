#! /bin/bash
sequence=(BasketballPass BQMall BQSquare ChinaSpeed KristenAndSara PartyScene SlideShow Traffic)
seqnum=${#sequence[@]}
source /etc/profile
for((i=0;i<seqnum;i++));do
bash /njulk/HEVC/allModels/sequences/mpmRunori.sh ${sequence[i]}
done
