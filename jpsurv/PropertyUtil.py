import ConfigParser

class PropertyUtil:
  cp = ConfigParser.SafeConfigParser()
  cp.optionxform=str
  cp.read(r"config.ini")
  properties = {}

  for section in cp.sections():
    for option in cp.options(section):
      properties[option]=cp.get(section,option)
  print properties
  @staticmethod
  def getAttribute(att):
    return PropertyUtil.properties[att]

  def getAsBoolean(self, path):
    return self.properties[path] in PropertyUtil.TRUTHY

  def getAsInt(self, path):
    return int(self.properties[path])

  def getAsFloat(self, path):
    return float(self.properties[path])

  def getAsString(self, path):
    return self.properties[path]


