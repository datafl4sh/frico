# FRICO examples

This folder contains some examples of the usage of FRICO, in particular:

* `biquad` In this folder you will find the model of a biquad antenna, which is quite popular in UHF/SHF amateur radio and between WiFi enthusiasts. It is the most complete example and you should start from it.
* `can` This is a variant of the infamous "Cantenna", or "Pringles antenna" and is sized according the dimensions found in the book "Le antenne riceventi e trasmittenti" by "Nuova Elettronica" (RIP).
* `dipole` A simple dipole model, mostly to check that the computed values are in the ballpark of the theoretical ones, or to compare with NEC.
* `yagi` A simple 5-element Yagi sized according to "Nuova Elettronica"'s book. Allows easy comparisons with NEC.

Before running the examples, please make sure to understand how FRICO outputs
data. This is described in the FRICO man page or in the HTML/PDF files included
under `doc/`.
