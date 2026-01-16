eps_arr=(0.1 0.05 0.01 0.005 0.001 0.0005)
k_arr=(5 10 15 20 25 30)



g++ -march=core2 -ffast-math -use_fast_math -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o lzpush SFMT.c main.cpp -g

mkdir result
mkdir result/dblp
cd result/dblp
cd ../..

./lzpush -d dblp -algo GEN_QUERY -n 20
./lzpush -d dblp -algo groundtruth -n 10
for ((i=0; i<${#eps_arr[@]}; i++))
do
    echo "./lzpush -d dblp -algo lanczos -n 10 -k ${k_arr[$i]}" |bash;
    echo "./lzpush -d dblp -algo lanczos_push -n 10 -k 20 -e ${eps_arr[$i]}" |bash;
    echo "./lzpush -d dblp -algo plot_results -n 10" |bash;
done