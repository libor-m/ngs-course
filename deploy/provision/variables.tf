variable "base_name" {
  description = "Base name for all Azure resources."
  type        = string
  default     = "ngs-course-2026"
}

variable "location" {
  description = "Azure region for provisioning."
  type        = string
  default     = "austriaeast"
}

variable "admin_username" {
  description = "Admin username for the VM."
  type        = string
  default     = "liborm"
}

variable "ssh_public_key_path" {
  description = "Path to the SSH public key used for VM access."
  type        = string
}
