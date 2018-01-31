from bsc.LocalDatabase import LocalDatabase

class LocalDatabaseDecorator(LocalDatabase):

    def createTable(self, tableName, tableDetail):
        
        # Create the table
        command = "CREATE TABLE %s(%s)"

        pass

    def deleteTable(self):
        pass

if __name__ == "__main__":

    pass