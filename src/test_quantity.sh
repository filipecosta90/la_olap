#!/bin/bash

~/postgres/bin/psql -U a57816 postgres -f /home/a57816/CPD_PI_OLAP/src/quantity.sql >> "/home/a57816/CPD_PI_OLAP/src/quantity_pgres.tbl"
