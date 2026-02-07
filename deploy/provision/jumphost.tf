locals {
  jumphost_name_prefix = "${var.base_name}-jumphost"
  jumphost_nsg_name    = "${local.jumphost_name_prefix}-nsg"
  jumphost_pip_name    = "${local.jumphost_name_prefix}-pip"
  jumphost_nic_name    = "${local.jumphost_name_prefix}-nic"
  jumphost_vm_name     = "${local.jumphost_name_prefix}-vm"
  jumphost_osdisk_name = "${local.jumphost_name_prefix}-osdisk"
}

resource "azurerm_network_security_group" "jumphost" {
  name                = local.jumphost_nsg_name
  location            = azurerm_resource_group.main.location
  resource_group_name = azurerm_resource_group.main.name

  security_rule {
    name                       = "${local.jumphost_name_prefix}-ssh"
    priority                   = 1001
    direction                  = "Inbound"
    access                     = "Allow"
    protocol                   = "Tcp"
    source_port_range          = "*"
    destination_port_range     = "22"
    source_address_prefix      = "*"
    destination_address_prefix = "*"
  }

  security_rule {
    name                       = "${local.jumphost_name_prefix}-https"
    priority                   = 1002
    direction                  = "Inbound"
    access                     = "Allow"
    protocol                   = "Tcp"
    source_port_range          = "*"
    destination_port_range     = "443"
    source_address_prefix      = "*"
    destination_address_prefix = "*"
  }
}

resource "azurerm_public_ip" "jumphost" {
  name                = local.jumphost_pip_name
  location            = azurerm_resource_group.main.location
  resource_group_name = azurerm_resource_group.main.name
  allocation_method   = "Static"
  sku                 = "Standard"
}

resource "azurerm_network_interface" "jumphost" {
  name                = local.jumphost_nic_name
  location            = azurerm_resource_group.main.location
  resource_group_name = azurerm_resource_group.main.name

  ip_configuration {
    name                          = "${local.jumphost_name_prefix}-ipcfg"
    subnet_id                     = azurerm_subnet.main.id
    private_ip_address_allocation = "Dynamic"
    public_ip_address_id          = azurerm_public_ip.jumphost.id
  }
}

resource "azurerm_network_interface_security_group_association" "jumphost" {
  network_interface_id      = azurerm_network_interface.jumphost.id
  network_security_group_id = azurerm_network_security_group.jumphost.id
}

resource "azurerm_linux_virtual_machine" "jumphost" {
  name                            = local.jumphost_vm_name
  resource_group_name             = azurerm_resource_group.main.name
  location                        = azurerm_resource_group.main.location
  size                            = "Standard_B2ts_v2"
  admin_username                  = var.admin_username
  disable_password_authentication = true

  network_interface_ids = [
    azurerm_network_interface.jumphost.id,
  ]

  admin_ssh_key {
    username   = var.admin_username
    public_key = file(var.ssh_public_key_path)
  }

  os_disk {
    name                 = local.jumphost_osdisk_name
    caching              = "ReadWrite"
    storage_account_type = "Premium_LRS"
    disk_size_gb         = 30
  }

  source_image_reference {
    publisher = "Debian"
    offer     = "debian-13"
    sku       = "13"
    version   = "latest"
  }
}
