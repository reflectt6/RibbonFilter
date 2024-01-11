//
// Created by Rain Night on 2024/1/10.
//

#include "org_example_JNIAdvancedDemo.h"

/*
 * Class:     org_example_JNIAdvancedDemo
 * Method:    nativeAdd
 * Signature: (Lorg/example/NameReceiver;Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_org_example_JNIAdvancedDemo_nativeAdd
        (JNIEnv * env, jclass jc, jobject jo, jstring js) {
    // 获取类引用
    jclass nrClass = env->GetObjectClass(jo);

    // 获取字段ID
    jfieldID maxCntField = env->GetFieldID(nrClass, "maxCnt", "I");
    jfieldID countField = env->GetFieldID(nrClass, "count", "I");

    // 获取字段值
    jint maxCnt = env->GetIntField(jo, maxCntField);
    jint count = env->GetIntField(jo, countField);

    // 在本地代码中使用字段值
    printf("Value of maxCnt: %d\n", maxCnt);
    printf("Value of count: %d\n", count);

    // 调用Java方法
    jmethodID addNameMethod = env->GetMethodID(nrClass, "addName", "(Ljava/lang/String;)V");
    env->CallVoidMethod(jo, addNameMethod, js);
}