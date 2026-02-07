#!/usr/bin/env bash
set -euo pipefail

# Required tools:
# - Azure CLI: az
# - Terraform: terraform
# - OpenSSH tools: ssh-keygen, ssh
#
# Usage (from repo root):
#   ./deploy/provision.sh
#
# After provisioning, SSH into the VM with:
#   ssh -i deploy/secrets/id_ed25519 <admin_username>@<public_ip>
#
# The public IP is printed at the end and available via:
#   terraform -chdir=deploy/provision output -raw public_ip

PROVISION_DIR="deploy/provision"
SECRETS_DIR="deploy/secrets"

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required tool: $1" >&2
    exit 1
  fi
}

require_cmd terraform


mkdir -p "$SECRETS_DIR"

# create azure principal
az ad sp create-for-rbac `
    --name "sp-ngs-course-2026-tf" `
    --role Contributor `
    --scopes /subscriptions/00f00dbd-7977-4267-9173-1984cc696cfb `
    --sdk-auth |
 jq -r '
"export ARM_CLIENT_ID=\"\(.clientId)\"
export ARM_CLIENT_SECRET=\"\(.clientSecret)\"
export ARM_TENANT_ID=\"\(.tenantId)\"
export ARM_SUBSCRIPTION_ID=\"\(.subscriptionId)\""' \
> $SECRETS_DIR/auth-azure.sh

source $SECRETS_DIR/auth-azure.sh

cp ~/.ssh/id_ed25519.pub "$SECRETS_DIR/id_ed25519.pub"

read -r -p "Admin username for the VM: " ADMIN_USERNAME

ADMIN_USERNAME=liborm

pushd "$PROVISION_DIR" >/dev/null

terraform init

# beware, the file paths are different here
terraform plan \
  -var "admin_username=$ADMIN_USERNAME" \
  -var "ssh_public_key_path=../secrets/id_ed25519.pub" \
  -out tfplan

terraform apply "tfplan"

PUBLIC_IP="$(terraform output -raw public_ip)"
JUMP_IP="$(terraform output -raw jumphost_ip)"

popd >/dev/null

ssh $ADMIN_USERNAME@$PUBLIC_IP

# TODO: add DNS A record
# and then
ssh liborm@ngs-course.behavio.dev
