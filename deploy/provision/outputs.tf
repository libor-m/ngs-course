output "public_ip" {
  description = "Public IP address of the VM."
  value       = azurerm_public_ip.main.ip_address
}

output "admin_username" {
  description = "Admin username used for SSH."
  value       = var.admin_username
}

output "resource_group" {
  description = "Resource group name."
  value       = azurerm_resource_group.main.name
}
