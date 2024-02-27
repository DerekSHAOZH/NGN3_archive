control="day0"
condition="day3"
sample=$condition"_background_"$control
data="/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMATS/data/"
gtf="/project/PausingDynamics/GeneralResources/Genecode29/gencode.v29.annotation.gtf"
out="/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMATS/Results/"$sample"/"
tmp="/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMATS/tmp/"$sample"/"
rmats="/home/gajos/Programs/rMATS/rmats.py"
len=100

python3 $rmats --b1 $data$control".txt" --b2 $data$condition".txt" --gtf $gtf -t paired --readLength $len --nthread 30 --od $out --tmp $tmp
