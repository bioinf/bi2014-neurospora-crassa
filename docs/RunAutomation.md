
Автоматизация запуска программ
====================================
Пусть имеется набор утилит (_BLAST+_, _Exonerate_, _GeneWise_, что угодно),
которые мы умеем запускать "ручками" и получать требуемый результат
(поиск выравнивания, сборка и прочее). Здесь "ручками" подразумевает запуск нужных программ
в консоли, передача им необходимых аргументов и обработка их (программ) вывода.

Зачем может понадобиться автоматизация?
----------------------------------------------
Рассмотрим, какие бонусы нам может дать автоматизация.

 *  __Контроль параметров запуска__
    Предположим, что программа для выравнивания имеет некоторые аргументы,
    регулирующие ее выполнение (кроме, разумеется, непосредственно входной и референсной
	последовательности). Это может быть алгоритм выравнивания, таблица _score_-ов, число
	потоков исполнения, порог ошибки, и т. п.

    Чтобы оценить масштаб "бедствия" посмотрим на число возможных аргументов _Exonerate_.

    ~~~bash
    $ exonerate --help
    ~~~

    Автоматизация позволит каждый последующий раз запускать программу ровно с теми же
	(или слегка измененными параметрами). Нет необходимости подбирать нужные параметры снова или
	переписывать их с блокнота.

    При этом мы можем легко изменить эти параметры или сделать несколько запусков с неким шагом
	какого-либо из параметров.
	
 *  __Работа с большим набором данных__
    Пусть имеется гигантский набор входных данных. Для каждого элемента набора
	необходимо запустить некоторые программы и, возможно, как-то дополнительно обработать результат.
	Ясно, что с ростом числа элементов в наборе, выполнить задачу "ручками" будет невозможным.

    Например, имеется $10^3$ различных белков (аминокислотных последовательностей), и требуется
	найти соответствующие им гены в некоторой референсной последовательности нуклеотидов. Если
	автоматизировать запуск программы с помощью, например, _python_, то для выполнения
	задачи достаточно лишь цикла `for`.

 *  __Создание конвеера из программ__
    Иногда для получения нужного результата недостаточно запустить одну утилиту и прочитать ее вывод.
	Может потребоваться взаимодействие нескольких утилит, например:
	 1.  Нужно запустить __a__, получить вывод `a.out`
	 2.  Запустить для `a.out` __b__ и __c__, получить `b.out`, `c.out`.
	 3.  Используя еще один __d__, который требует `c.out` и `b.out`, получить уже финальный
         результат.

    Заметим, что если запуск всех программ был бы выражен через _python_, то
	для выполнения задачи останется:

	~~~python
    def complex_procedure(input):
	    a = do_a(input)
		b = do_b(a)
		c = do_c(a)
		return do_d(b, c)
    ~~~

    А если подобную процедуру необходимо повторить много раз, то без автоматизации не обойтись.

Как взаимодействовать с программой?
-----------------------------------------
Для начала, чтобы автоматизировать запуск программ, необходимо четко представлять, что
этот запуск собой представляет. Все взаимодействие с программой можно условно разделить на четыре
части:

 *  __Аргументы командной строки__ --- это конкретные параметры запуска:

    ~~~bash
    $ find ~/Downloads -iname "the_answer_to_life*"
    ~~~

    Здесь `.`, `-iname` и  `"the_answer_to_life*"` суть аргументы командной строки.

	
 *  __Вход программы__ --- это, к примеру, пользовательский ввод или перенаправление из файла
    с использованием соответствующих возможностей _Bash_.
	
    ~~~bash
    $ sort < data.txt
	$ sudo su
	Password:
    ~~~

    Содержимое `data.txt` является входом программы.
   
 *  __Выход программы__ --- это все то, что программа печатает на консоль в процессе своей работы.

    ~~~bash
	$ which less
	/usr/bin/less
    ~~~
    
    Выходом программы является строка `/usr/bin/less`. 

 *  __Файлы__ --- это все те файлы, в которые программа пишет или которые читает.

    ~~~bash
	$ genblasta -q sample.proteins -s neurospora.fasta -o genblasta.out
	~~~

    В данном примере речь идет о файлах `sample.proteins`, `neurospora.fasta`, `genblasta.out`.

Разобравшись с тем, в чем именно заключается взаимодействие пользователя с программой,
перейдем к реализации этого взаимодействия с помощью _python_.

Работа с файлами в _python_
------------------------------
Итак, представим текстовый редактор (например _Gedit_). Что самое простое умеет делать каждый текстовый редактор?

__Открывать и закрывать файлы.__

Как ни странно, _python_ умеет делать тоже самое. При этом точно
также, как в текстовом редакторе каждому открытому файлу соответствует свое окошко (вкладка, буфер),
в _python_ каждому открытому файлу соответствует специальный идентификатор, который создается
функций `open()`:

~~~python
f = open('data.txt', 'r')
f.close()
~~~

Функция `open()` принимает два аргумента-строки: __путь к файлу__ и __режим работы__.
Существует четыре режима работы:

 *  __`'r'`__ --- только чтение;
 *  __`'w'`__ --- только запись;
 *  __`'r+'`__ --- и чтение, и запись;
 *  __`'a'`__ --- запись в конец файла (оставляя предыдущее содержимое).

__Внимание!__ Файл необходимо всегда закрывать (`close()`) после работы с ним, иначе
можно столкнуться с ситуацией, когда тот же самый файл окажется открыт еще раз нами же или
в другом "текстовом редакторе".
A такая ситуация приводит к проблемам несогласованных изменений даже в случае
с обычными текстовыми редакторами,  _python_ --- не исключение.

Что еще можно ожидать от любого текстового редактора? Возможность ввода новой информации
в и чтения имеющейся. Таким образом, если файловый идентификатор хранится в переменной `f`,
то можно совершать следующие базовые операции:

 *  `f.read()` --- возвращает строку со всем содержимым файла.
 *  `f.write(str)` --- дописывает содержимое строчки `str` в файл.
 *  `f.readline()` --- возвращает следующую (после последней считанной)
    строку файла вместе с символом окончания строки `'\n'`
    или пустую строку, если больше строк нет.

Эти операции читают/пишут последовательно, т. е., например, если в файле
`neurospora_crassa_or74a_12_supercontigs.fasta` две строчки:

~~~
>Supercontig_12.1 of Neurospora crassa OR74A
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
~~~

Тогда работа с файлом в интерпретаторе может выглядеть так:

~~~python
$ python
Python 3.3.5 (default, Aug 24 2014, 22:55:44) 
[GCC 4.7.3] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> fasta = open('test.fasta', 'r')
>>> fasta.readline()
'>Supercontig_12.1 of Neurospora crassa OR74A\n'
>>> fasta.readline()
'CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC\n'
>>> fasta.readline()
''
>>> fasta.close()
>>>
~~~

__"Финт ушами"__ Есть магия, с помощью которой можно избежать явного закрытия файла:

~~~python
with open(...) as f:
  # работаем с идентификатором файла f
~~~

При выходе из блока `with` файл будет автоматически закрыт.

Стандартные идентификаторы файлов
----------------------------------------

Обсудив работу с файлами, используем полученную информацию для
работы с вводом и выводом программы. Оказывается, можно считать,
что при вводе/выводе программа на самом деле читает/пишет из некоторых специальных файлов,
идентификаторы которых ей кто-то (операционная система, вызывающая программа) передает "свыше".
Назовем идентификатор, который соответствует входу программы, __stdin__,
а файл, который соответствует выводу программы, __stdout__.

Далее, если мы запускаем наш скрипт, например, так:

~~~bash
$ ./sweet.py
~~~

Тогда:

 * __stdin__ --- идентификатор, который соответствует вводу пользователя в консоль;
 * __stdout__ --- идентификатор, который соответствует экрану консоли.

Если же мы запускаем скрипт как:

~~~bash
$ ./sweet.py <inputfile >outputfile
~~~

То:

 *  __stdin__ --- идентификатор, соответствующий открытому в режиме "для чтения" файла `inputfile`;
 *  __stdout__ --- идентификатор, соответствующий открытому в режиме "для записи" файла `outputfile`.

Внутри _python_ доступ к этим файлам, например, спрятан в функциях `input()` и `print()`.
К этим файлам можно также получить прямой доступ через `sys.stdin` и `sys.stdout`:

~~~python
>>> import sys
>>> sys.stdin
<_io.TextIOWrapper name='<stdin>' mode='r' encoding='UTF-8'>
>>> sys.stdout.write("Hello\n")
Hello
6
>>> 
~~~

Не так важно знать про то, как получить доступ к этим спец-идентификаторам, важно понимать, что
ввод и вывод со стороны программы это ровно то же самое, что чтение/запись из файла.


Списки в _python_
-------------------
В _python_ есть тип данных _список_. Что же это такое?
Вспомним как выглядят списки в обычной жизни. Например, список покупок:

 1. Terra Mystica
 2. Power Grid
 3. Ticket to ride

Или еще список:

 1. fungus
 2. βακτήριον
 3. ἀρχαῖα
 4. virus
 5. animalia
 6. plantae

По сути, список ничего более чем последовательность произвольных элементов.
В _python_ список -- это ровно то же самое --- последовательность произвольных значений.
Синтаксис для создания списков следующий:

~~~python
games = ['Terra Mystica', 'Power Grid', 'Ticket to Ride']
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
heterogeneous = [1, 'rock', 2, 'paper', 3, 'scissors', 4, 'lizard', 5 'Spock']
~~~

Над списками можно производить следующие операции:

~~~python
>>> numbers[3] # получить элемент по номеру (начиная с нуля)
4
>>> numbers[1:3] # получить подсписок по двум индексам (второй индекс не включается)
[2, 3]
>>> numbers[1:3] = ['a', 'b', 'c'] # заменять подсписок на другой
>>> numbers
[1, 'a', 'b', 'c', 4, 5, 6, 7, 8, 9, 10]
>>> numbers[1:3] = [] # удалять элементы из списка
>>> numbers
[1, 'c', 4, 5, 6, 7, 8, 9, 10]
>>> [1, 2, 3] + [1, 2, 4] # объединять списки
[1, 2, 3, 1, 2, 4]
>>> len(numbers) # узнать длину списка
9
>>> for x in numbers:
...     # выполнить некоторые действия для всех элементов списка
...     print(x)
... 
~~~

Модуль `subprocess`
----------------------------
Вернемся к нашей задаче автоматизации запуска программ.
В самом начале были выделены четыре основных способа взаимодействия
с программой, при этом, как работать в _python_ с последним из них, мы уже обсудили.
Для работы с оставшимися тремя в _python_ имеется специальный модуль `subprocess`.

Нас интересует лишь функция `call()` этого модуля. Рассмотрим интересующие нас ее параметры:

~~~python
subprocess.call(args, stdin=None, stdout=None)
~~~

 *  __`args`__ Список строк, первым элементом которого является имя программы, которую мы хотим
    запустить, а последующие элементы --- аргументы командной строки для передачи данной
	команде.
 *  __`stdin`__ Необязательный (если отсутствует, используется наш же собственный ввод) именованный
    параметр, являющийся идентификатором открытого на чтение (как минимум) файла.
 *  __`stdout`__ Необязательный (если отсутствует, используется наш же собственный вывод) именованный
    параметр, являющийся идентификатором открытого на запись (как минимум) файла.

Возвращает функция `call()` код возврата (целое число) выполненной программы.
Как общее правило: нулевой код возврата соответствует успешному завершению, любой другой ---
какой-либо ошибке.

Пример:

~~~python
>>> import subprocess
>>> f = open('data.txt', 'r')
>>> return_code = subprocess.call(['sort', '-r'], stdin = f)
...
>>> return_code
0
>>> f.close()
~~~

Здесь снова используется программа `sort`, которая сортирует свой вход построчно и
выводит результат (опция `-r` включает обратный порядок сортировки).
После выполнения представленного кода содержимое файла `data.txt`
будет отсортировано в обратном порядке и выведено на консоль.

Биоинформатичный пример, запускающий `exonerate` и записывающий вывод в файл `exonerate.out`:

~~~python
def run_exonerate(proteins, sequence):
	with open('exonerate.out', 'w') as outfile:
		subprocess.call(['exonerate', '--model', 'protein2genome', proteins, sequence],
	                    stdout=outfile)
~~~

Переменные окружения
-------------------------
В качестве небольшого дополнения на уровне рецепта рассмотрим еще один параметр к `subprocess.call()`.

Наченем с того, что есть как минимум еще один неупомянутый
способ взаимодействия с программой --- переменные окружения, к которым программа также имеет доступ.
Для того, чтобы указать при запуске программы, какие же у нее в окружении будут
переменные, `subprocess.call` имеет именованный параметр __`env`__, который принимает т. н.
_словарь_ (пара ключ-значение). Рассмотрим пример работы со словарями (все вполне интуитивно):

~~~python
>>> marks = { 'John' : 'A', 'Peter' : 'B', 'Lisa' : 'F' }
>>> marks['Lisa'] = 'A'
>>> marks
{'Lisa': 'A', 'John': 'A', 'Peter': 'B'}
>>> marks['Peter']
'B'
>>> del marks['Peter']
>>> marks
{'Lisa': 'A', 'John': 'A'}
~~~

Для того, чтобы получить словарь соответствующий текущему окружению, нужно взять поле
`environ` из модуля `os`:

~~~python
import os
my_env = os.environ
~~~

Модификация окружения производится с помощью редактирования словаря:

~~~python
new_env = os.environ.copy()
new_env['PATH'] = '/home/me/soft/bin/:' + new_env['PATH']
# Здесь используется функция copy(), создающая копию словаря, т. к.
# иначе оказалось бы измененным наше собственное окружение
~~~

Передача нового окружения в `subprocess.call()` уже не представляет никаких сложностей:

~~~python
import subprocess
subprocess.call(['printenv'], env=new_env)
~~~

В результате последнего вызова, в консоль будет выведен список всех переменных окружения,
при этом значение `PATH` должно теперь содержать `/home/me/soft/bin/`.

Пример c запуском _GenBlastA_
-------------------------------
Известно, что _GenBlastA_ требует для своего выполнения старый _BLAST_.
Однако в нашем `PATH` лежит только путь к новому _BLAST+_. Вспоминая, что
директории в `PATH` проверяются слева направо, можно прийти к следующему способу запуска
_GenBlastA_: на время запуска дописать к _PATH_ слева путь к папке `bin` старого _BLAST_.

В _Bash_ это можно сделать, используя команду `env`, которая предназначена именно для
модификации окружения запускаемой программы:

~~~bash
$ env HI=ha printenv HI
ha
~~~

А в _python_ можно использовать разобранный способ:

~~~python
old_blast_env = os.environ.copy()
old_blast_env['PATH'] =  '/home/me/soft/oldblast/bin:' + old_blast_env['PATH']
subprocess.call(['genblasta', '-q', proteins, '-t', sequence, '-o', 'genblasta.out'],
                env=old_blast_env)
~~~
