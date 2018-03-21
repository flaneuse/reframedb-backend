# Flask JWT Auth

[![Build Status](https://travis-ci.org/realpython/flask-jwt-auth.svg?branch=master)](https://travis-ci.org/realpython/flask-jwt-auth)

## Want to learn how to build this project?

Check out the [blog post](https://realpython.com/blog/python/token-based-authentication-with-flask/).

## Want to use this project?

### Basics

1. Fork/Clone
1. Activate a virtualenv
1. Install the requirements

### Set Environment Variables

Update *project/server/config.py*, and then run:

```sh
$ export APP_SETTINGS="project.server.config.DevelopmentConfig"
```

or

```sh
$ export APP_SETTINGS="project.server.config.ProductionConfig"
```

Set a SECRET_KEY:

```sh
$ export SECRET_KEY="change_me"
```

Set the secret key for your Google ReCaptcha:

```sh
$ export RECAPTCHA_SECRET_KEY="change_me"
```

Set the secret key for your Google ReCaptcha:

```sh
$ export DATA_DIR="directory_with_data_files"
```

Set the mail username for user confirmation:

```sh
$ export API_MAIL_USERNAME="username"
```

Set the mail password:

```sh
$ export API_MAIL_PASSWORD="password"
```

### Create DB

Create the databases in `psql`:

```sh
$ psql
# create database repurpos_db;
# create database repurpos_db_test;
# \q
```

Create the tables and run the migrations:

```sh
$ python manage.py create_db
$ python manage.py db init
$ python manage.py db migrate
```

### Run the Application

```sh
$ python manage.py runserver
```

Access the application at the address [http://localhost:5000/](http://localhost:5000/)

> Want to specify a different port?

> ```sh
> $ python manage.py runserver -h 0.0.0.0 -p 8080
> ```

### Testing

Without coverage:

```sh
$ python manage.py test
```

With coverage:

```sh
$ python manage.py cov
```
