rm *.ans
scp -P 5522 -i ~/key/LAB_RSA.txt oba@140.123.104.79:/home/oba/Fidelity-Aware-Swapping-Ton-3rd/data/ans/data.zip .
unzip data.zip