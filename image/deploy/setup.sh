
# prepare passwords for users
# aptitude search '?provides(wordlist)'
sudo apt install wamerican

# generate some funny passwords
</usr/share/dict/words egrep "^[a-z]{5,8}$" |
  sort -R |
  paste -d' ' - - - |
  head -30 |
  nl -w2 -n'rz' |
  sed 's/^/user/' \
> image/deploy/secrets/users.tsv

export ANSIBLE_CONFIG="image/deploy/setup/ansible.cfg"

ansible-playbook \
  -i "image/deploy/setup/hosts.ini" \
  "image/deploy/setup/site.yml"

# set up the jumphost only
ansible-playbook \
  -i "image/deploy/setup/hosts.ini" \
  "image/deploy/setup/jumphost.yml"

# finish things that are at the end and replaying the whole book
# takes so looong
ansible-playbook \
  -i "image/deploy/setup/hosts.ini" \
  "image/deploy/setup/finish.yml"