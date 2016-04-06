#!/usr/bin/python

import sys
import sqlalchemy
from SQLiteConnection import engine, Session
from ModelClasses import *



filename = 'student_data.txt'

data = open(filename)

session = Session()
#student = Student()
for line in data:
	print(line)
	if line[0] == '#':
		continue
	else:
		#for the column names		
		
		line = line.strip('\n')
		line = line.split('|')        
			
		student = Student()
		student.first_name = line[0]
		student.last_name = line[1]
		session.add(student)
	
	for supervisor_room in line[3].split(', '):
		if supervisor_room != []		
		#supervisor,room = supervisor_room.split('/')

		try:
			one_supervisor = session.query(Supervisor).filter(Supervisor.last_name==supervisor_room).one()
		except sqlalchemy.orm.exc.NoResultFound:
			one_supervisor = Supervisor()
			one_supervisor.last_name = supervisor_room
			session.add(one_supervisor)
		except sqlalchemy.orm.exc.MultipleResultsFound:
			print ('There is more than one Doctor!')
			sys.exit(1)
		student.supervisors.append(one_supervisor)

	

#	try:
#		one_city = session.query(City).filter(City.last_name==line[2]).one()
#	except sqlalchemy.orm.exc.NoResultFound:
#		one_city = City()
#		one_city.last_name = line[2]
#		session.add(one_city)
#	except sqlalchemy.orm.exc.MultipleResultsFound:
#		print "There is more than one city!"
#		sys.exit(1)

#	student.cities.append(one_supervisor)	



session.commit()

engine.dispose() # cleanly disconnect from the database
sys.exit(0)
