for i in $(seq 73270 1 73317); do 
   echo ${i}
   scancel ${i} 
done

