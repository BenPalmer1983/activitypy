from activity import activity

class end:

  def run(code):
    if(code == 'invalid_projectile'):
      end.error()
      print(
"""Invalid projectile code selected.
Accepted codes:
p = proton
d = Deuteron
t = Triton
he3 = Helium-3
a = alpha particle
""")
    else:
      print("Unknown error")
    exit()

  def error():
     print(
"""
################################
Program terminated with an error
Time: """ + str(activity.time()) + """
################################

""")
