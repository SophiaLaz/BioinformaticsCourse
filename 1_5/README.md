## Домашнее задание по пятой лекции.
### Практическая часть
По написанию кода и использованию библиотек все как в 
первом ДЗ. Для тестов можно использовать 
последовательности, результат выполнения алгоритма на 
которых можно проверить вручную.
1. **Жадный алгоритм множественного выравнивания** (7 баллов) 
Реализуйте алгоритм, который принимал бы на вход 
массив строк, штраф за удаления, вставки и несовпадения,
а так-же цену совпадений. А возвращал бы множественное 
выравнивание. На первом шаге алгоритм должен выбрать 
две самые близкие по расстоянию Левенштейна строки и 
заменить их консенснусной строкой. При следующих шагах 
алгоритма выравниваться между собой могут так же и 
консенснусные строки. При этом стоит хранить для каждой 
строки не только ее саму но и профиль множественного 
выравнивания, чтобы в итоге правильно пересчитывать 
консенсус.  
Результат работы алгоритма - массив строк, 
соответствующий некоторому множественному выравниванию.