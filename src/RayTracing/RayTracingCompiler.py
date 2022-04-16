import re

def Compiler(Str):
  Str = re.sub(r'\n\s*\n', r'\n', Str)
  Str = re.sub(r'\n\s*', r'\n', Str)

  Str = Str.split("\n")
  Codes = "void Model(RayTracing& ray) {\n\tMaterial* material;\n\t"

  for i in range(len(Str)):
    s = Str[i].split("\t")

    if(s[0] == '' or s[0][0] == '/'):
      continue

    # Material
    if(s[0] == 'Material'):
      Codes = Codes + "material = new Material"

      for j in range(1, len(s)):
        Codes = Codes + "; material->" + s[j]

      Codes = Codes + ";\n\t"
      continue

    Codes = Codes + "ray.objTree.add" + s[0] + '('
    Codes = Codes + s[1] + ', material'
    if (len(s) == 3):
      Codes = Codes + s[2]
    Codes = Codes + ");\n\t"

  Codes = Codes + "\n}\n"
  return Codes



if __name__ == '__main__':
    fileName = "D:/te.txt"
    file = open(fileName,"r", encoding='utf-8')
    Str = file.read()
    file.close()

    Str = Compiler(Str)

    file = open("D:/sa.c","w+", encoding='utf-8')
    file.write(Str)
    file.close()