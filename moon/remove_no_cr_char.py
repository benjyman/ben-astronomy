filename = 'earthshine_fmlist.py'
with open(filename, 'rb+') as f:
    content = f.read()
    f.seek(0)
    f.write(content.replace(b'\r', b''))
    f.truncate()
