BayesRadonReconst
=========================================

Bayes推定に基づく Radon 変換の再構成プログラム．
R によって記述子ているが，Radon 変換のシミュレータと
フィルタ逆投影に関しては C 言語で記述している．

実験方法
------------------------

1. まず， Cobj において `make` を実行し，Radon 変換のシミュレータと
   フィルタ逆投影プログラムをコンパイルする．

2. MakeData.R を実行することで，サイノグラムデータが，
   `P1Reconst.Rdata` にセーブされる
	   
3. RecDemo.R を実行すると再構成のデモが行われる．
   描画される画像は，通常の微分フィルタを用いた FBP
   HAN窓をフィルタとした FBP
   MRF を事前分布とし，ハイパーパラメータ推定まで行ったフィルタで再構成させた Bayes FBP の画像が描画される
   

<img width="1434" alt="recdemo" src="https://cloud.githubusercontent.com/assets/5178986/24691284/8bb4c31e-1a0d-11e7-9a55-8828a77f37f5.png">
