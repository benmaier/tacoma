#!/bin/bash

convert new_logo.png -set colorspace Gray -separate -average new_logo_grey.png
convert new_logo_small.png -set colorspace Gray -separate -average new_logo_small_grey.png
