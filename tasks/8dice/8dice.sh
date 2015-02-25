# sample problem
http://leancrew.com/all-this/2015/02/counting-heads/
Eight people are sitting around a circular table. Each has a coin. They all flip their coins. What is the probability that no two adjacent people will get heads?

(echo "obase=2;"; printf "%d\n" {0..255}) | bc |  # generate numbers 0-255 in binary
  sed 's/^/0000000/' | rev | cut -c1-8 | rev |    # pad the output to 8 characters
  sed p | paste - - | tr -d "\t" |                # double the pattern on each line to simulate ring
  grep -v 11 |                                    # ignore all patterns where two heads (1) are next to each other
  wc -l                                           # count the rest

# modular variant
generate () { (echo "obase=2;"; echo {0..255} | tr " " "\n") | bc ;}
generate () { (echo "obase=2;"; echo {0..255} | xargs -n1) | bc ;}
generate () { (echo "obase=2;"; printf "%d\n" {0..255}) | bc ;}

pad () { awk '{print substr(sprintf("0000000%s", $0), length);}' ;}
pad () { sed 's/^/0000000/' | rev | cut -c-8 | rev ;}
pad () { sed 's/^/0000000/' | egrep -o '.{8}$' ;}

ring () { sed p | paste - - | tr -d "\t" ;}
ring () { sed 'h;G;s/\n//' ;}

generate | pad | ring | grep -v 11 | wc -l

# final version
(echo "obase=2;"; printf "%d\n" {0..255}) | bc |  # generate numbers 0-255 in binary
  sed 's/^/0000000/' | egrep -o '.{8}$' |         # pad the output to 8 characters
  sed 'h;G;s/\n//' |                              # double the pattern on each line to simulate ring
  grep -v 11 |                                    # ignore all patterns where two heads (1) are next to each other
  wc -l                                           # count the rest
