#!/bin/sh

for p in $(ls ../problems/*/example.bc | sed 's/\/example\.bc//g')
do
	for ex in $p/example*.bc; do
        echo '===================================='
		file=`basename $ex .bc`
        cmd="$p/biqcrunch $ex $p/biq_crunch.param"
        echo $cmd
        $cmd | grep 'value' > $p/${file}.out
        if [ -z "`diff -q $p/${file}.out $p/${file}.ans`" ]
        then echo `cat $p/${file}.out` ... OK
        else echo `cat $p/${file}.out` ... FAIL
        fi
        rm -f $p/${file}.bc.output* $p/${file}.out
	done
done
echo '===================================='

