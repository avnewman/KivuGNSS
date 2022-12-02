#!/usr/bin/env bash

message="Auto-update: $(date)"
git add . 
git commit -m "$message"
git push


