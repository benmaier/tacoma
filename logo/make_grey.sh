#!/bin/bash

convert new_logo.png -set colorspace Gray -separate -average new_logo_grey.png
