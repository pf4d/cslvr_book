#!/bin/bash

for f in $(find . -name '*.png'); do
  fn=${f::-4};
  convert $fn.png -quality 70 -sharpen 0x1.0 $fn.jpg;
  rm $fn.png;
done



