# Molecular Project

## Overview

The Molecular Project is designed to manage and manipulate molecular structure data stored in XYZ files using a SQLite
database. This project includes functionality to create a database, insert data from XYZ files, and view the stored
data. Additionally, it provides a RESTful API to interact with the molecular data and includes unit tests to ensure the
correctness of database operations and API endpoints.

## Features

### Database Management:

Create and manage a SQLite database for storing molecular data.

### Data Insertion:

Insert molecular structure data from XYZ files into the database.

### Data Retrieval:

View and query the stored molecular data.

### API Endpoints:

Access and manipulate molecular data through a RESTful API.

### Unit Testing:

Ensure the accuracy and reliability of database operations and API endpoints.

## Requirements

Python 3.7+  
FastAPI  
SQLite  
Requests (for testing)

## Installation

### Clone the repository:

```shell
git clone https://github.com/yourusername/Molecular_Proj.git
cd Molecular_Proj
```

### Create a virtual environment and activate it:

```shell
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

### Install the dependencies:

```shell
pip install -r requirements.txt
```

### Set up your database and directories in settings.py:

```shell
db_path = "path/to/your/database.db"
directory_path = "path/to/your/xyz/files"
```

## Running the Application

### To start the FastAPI application, run:

```shell
python main.py
```

### The API will be available at

The API documentation is available at [Molecular_project_API_docs](http://localhost:8000/docs).

## Running the Tests

To run the tests, use:

```shell
pytest api_tests.py
```

This will execute the unit tests for the API endpoints and ensure that all functionalities are working correctly.