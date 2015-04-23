Starter session
===============

This session will give you all the basics that you need 
to smoothly move around when using a UNIX system (in the text mode!).

**Use multiple windows (and be safe when the network fails)**

First, type ``screen`` in your terminal::

  screen

Screen creates the first window for you. To create another one press 
``ctrl+a c``. To switch between the windows press ``ctrl+a space``.

**Check what the computer is doing**

Run ``htop`` in one of your screen windows::

  htop

Htop displays CPU and memory utilization of the (virtual) computer. Continue your 
work in another window (``ctrl+a space``). You can switch back to the htop window to
monitor progress of some lengthy operation.
