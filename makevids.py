def makevid_onecase(folderprepath, case, field):
    import subprocess
    # case = "exp1_coarse3"; figprepath = f"/moonbow/ascheb/escape/idealized/{case}/figures/"
    # field = "thetaxc"
    vidscript = open("make_vid_gif_commented.sh", "r+")
    lines = []
    for line in vidscript:
        if line.startswith("FOLDER_NAME"):
            # print("Found a folder!")
            line = f'FOLDER_NAME={folderprepath}{case}/figures/{field}\n'
        lines.append(line)
        lines.append("")
    
    vidscript.close();
    # print(lines)
    # raise Exception("Got the right lines; now need to write and run")
    newvidscript = open(f"make_vid_gif_commented.sh", "w")
    newvidscript.writelines(lines)
    newvidscript.close()
    print("Running!")
    subprocess.run(f"make_vid_gif_commented.sh")
    return None

def makevidcomb(folderprepath, field):
    import subprocess
    # case = "exp1_coarse3"; figprepath = f"/moonbow/ascheb/escape/idealized/{case}/figures/"
    # field = "thetaxc"
    vidscript = open("make_vid_gif_commented.sh", "r+")
    lines = []
    for line in vidscript:
        if line.startswith("FOLDER_NAME"):
            # print("Found a folder!")
            line = f'FOLDER_NAME={folderprepath}/{field}\n'
        lines.append(line)
        lines.append("")
    
    vidscript.close();
    # print(lines)
    # raise Exception("Got the right lines; now need to write and run")
    newvidscript = open(f"make_vid_gif_commented.sh", "w")
    newvidscript.writelines(lines)
    newvidscript.close()
    print("Running!")
    subprocess.run(f"make_vid_gif_commented.sh")
    return None



