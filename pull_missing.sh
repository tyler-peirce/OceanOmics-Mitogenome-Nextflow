#!/bin/bash

# Define source and destination directories
SOURCE_DIR="../_mounted"
DEST_DIR="../missing"

# Ensure the destination directory exists
mkdir -p "$DEST_DIR"

# Read each directory name from missing_list.txt and copy it
while IFS= read -r dir_name; do
    echo "cp -r "$SOURCE_DIR/$dir_name" "$DEST_DIR/""
    # Check if the directory exists in the source directory
    if [ -d "$SOURCE_DIR/$dir_name" ]; then
        echo "Copying $dir_name..."
        rclone copy "$SOURCE_DIR/$dir_name" "$DEST_DIR/$dir_name" --include "*.gz" --checksum --progress
        
    else
        echo "Warning: $dir_name not found in $SOURCE_DIR"
    fi
done < missing_list.txt

echo "Copy process completed!"
