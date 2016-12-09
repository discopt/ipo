class NonConstError:
    def init(self, value):
        self.value = 'This is no const value: '+value

    def __str__(self):
        return repr(self.value)
