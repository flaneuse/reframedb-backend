# project/server/config.py

import os
basedir = os.path.abspath(os.path.dirname(__file__))
postgres_local_base = 'postgresql://laurahughes:@localhost/'
database_name = 'repurpos_db'


class BaseConfig:
    """Base configuration."""
    SECRET_KEY = os.getenv('SECRET_KEY', 'my_precious')
    SECRET_PASSWORD_SALT = os.getenv('SECRET_PASSWORD_SALT', 'my_precious_two')
    DEBUG = False
    BCRYPT_LOG_ROUNDS = 13
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # Mail settings

    MAIL_SERVER = os.getenv('SES_MAIL_SERVER', 'email-smtp.us-east-1.amazonaws.com' )
    MAIL_PORT = 587
    MAIL_USE_TLS = True
    MAIL_USE_SSL = False
    MAIL_USERNAME = os.getenv('SES_SMTP_USERNAME')
    MAIL_PASSWORD = os.getenv('SES_SMTP_PASSWORD')
    MAIL_DEFAULT_SENDER = 'help@reframedb.org'



class DevelopmentConfig(BaseConfig):
    """Development configuration."""
    DEBUG = True
    BCRYPT_LOG_ROUNDS = 4
    SQLALCHEMY_DATABASE_URI = postgres_local_base + database_name
    FRONTEND_URL = 'http://localhost:4200/'


class TestingConfig(BaseConfig):
    """Testing configuration."""
    DEBUG = True
    TESTING = True
    BCRYPT_LOG_ROUNDS = 4
    SQLALCHEMY_DATABASE_URI = postgres_local_base + database_name + '_test'
    PRESERVE_CONTEXT_ON_EXCEPTION = False
    FRONTEND_URL = 'http://localhost:4200/'


class ProductionConfig(BaseConfig):
    """Production configuration."""
    SECRET_KEY = 'my_precious'
    DEBUG = False
    SQLALCHEMY_DATABASE_URI = postgres_local_base + database_name
    FRONTEND_URL = 'https://reframedb.org/'
