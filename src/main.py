import uvicorn
from database import create_database, delete_database

if __name__ == "__main__":
    create_database()
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)

# delete_database()
