import os


# A function to create folders
def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print("Created a folder: " + path) 
    return path

# A function to create links
def create_link(original, target):
    os.system("ln -s " + original + " " + target)
    print("Created a link for " + original + " at " + target)
    return target

# A function to copy files
def copy_file(original, target):
    os.system("cp " + original + " " + target)
    print("Copied the file " + original + " to " + target)
    return target

# Change the value of a specific entry in a text file
def change_entry(flag, value, file):
    os.system("sed -i 's/" + flag + "/" + value + "/' " + file)
    print("Changed the value of the entry " + flag + " to " + value + " in " + file)
    
# change directory
def change_directory(dir):
    os.chdir(dir)
    print("Current working directory: " + os.getcwd())
    return dir
    
# A function to run a system command
def system_com(command):
    print("running command: " + command)
    os.system(command)